!***********************************************************************
!  wrapper.f90 — a thin, ccall-friendly shim around ACM Calgo 766
!  VECTOR_PADE (Cabay, Jones & Labahn 1997, ACM TOMS 23(1), 91-110,
!  DOI 10.1145/244768.244790).
!
!  WHY THIS FILE EXISTS
!  --------------------
!  The shipped Calgo-766 routine VECTOR_PADE has a Fortran-90 interface
!  built entirely from *assumed-shape* arrays, and it depends on a host
!  module `working_area_VECTOR_PADE` of allocatable scratch arrays that
!  the *caller* must provide.  An assumed-shape dummy argument is passed
!  as a compiler-private array descriptor, not a bare pointer, so it
!  cannot be `ccall`-ed from Julia in any stable way.
!
!  This file therefore supplies (1) the `working_area_VECTOR_PADE`
!  module VECTOR_PADE needs, and (2) a single `bind(C)` entry point
!  `vpade_type2` whose entire ABI is plain `real*8` / `int32` arrays
!  passed by reference in column-major order — exactly what Julia's
!  `ccall` speaks.  Inside, it allocates fixed-shape arrays, calls
!  VECTOR_PADE, applies the documented `gamma_star` row scaling, and
!  copies the type-II (simultaneous-Pade / shared-denominator) row-0
!  result back to the caller.
!
!  CALGO-766 TYPE-II CONVENTION (see vector_pade.f90 header, lines
!  ~30-50 and lines 135-170 of Src/Sp/vector_pade.f90)
!  ------------------------------------------------------------------
!  Given k+1 power series A(0..k) and a degree vector n(0..k) with
!  ||n|| = sum(n), VECTOR_PADE returns the scaled simultaneous-Pade
!  system S_star.  Row 0 of S_star is the type-II approximant:
!
!      S_star(:,0,beta) * A(0)  -  S_star(:,0,0) * A(beta)
!                                                  = O(z**(||n||+1))
!
!  for beta = 1..k (this is exactly `build_delta_T_star`, alpha=0).
!  Hence with A(0) := the constant series 1 and A(beta) := f_beta,
!
!      Q       := S_star(:,0,0)        (degree ||n|| - n(0))
!      P_beta  := S_star(:,0,beta)     (degree ||n|| - n(beta))
!      f_beta  ~=  P_beta / Q          (to order z**(||n||+1)).
!
!  S_star as returned is *scaled*; the normalised system divides row i
!  by gamma_star(i).  Row 0 / gamma_star(0) is therefore the genuine
!  shared-Q approximant.  The Calgo coefficient order is low-to-high
!  (z**0 first), the same as SharedPade's `(numerators, denominator)`.
!  Calgo does NOT impose Q(0)=1; the caller (capture.jl) renormalises.
!
!  ABI of vpade_type2 (all arguments by reference; arrays column-major)
!  ------------------------------------------------------------------
!    k      in   int32             number of components d (k+1 series)
!    nvec   in   int32(k+1)        degree vector n(0..k)
!    nrm    in   int32             ||n|| = sum(nvec); A has nrm+1 rows
!    Amat   in   real8(nrm+1,k+1)  series; col 0 = ref series A(0),
!                                  col j = f_j.  Column-major: element
!                                  (l,j) is the z**l coeff of A(j).
!    tau    in   real8             Calgo stability tolerance
!    Qout   out  real8(nrm+1)      shared denominator Q, low-to-high,
!                                  zero-padded to nrm+1 entries
!    Pout   out  real8(nrm+1,k)    numerators P_1..P_k, column j-1 is
!                                  P_j low-to-high, zero-padded
!    flag   out  int32             Calgo error flag (0 ok, 1 ill-cond,
!                                  2 singular, 3 bad input)
!***********************************************************************

      module working_area_VECTOR_PADE
!        Allocatable scratch arrays required by VECTOR_PADE itself.
         real, dimension(:,:,:), allocatable :: S_hat, New_S
         real, dimension(:,:),   allocatable :: T
         real, dimension(:,:,:), allocatable :: S_star_hat, New_S_star
         real, dimension(:,:,:), allocatable :: T_star
      end module working_area_VECTOR_PADE


      subroutine vpade_type2(k, nvec, nrm, Amat, tau, Qout, Pout, flag) &
                 bind(C, name="vpade_type2")

      use, intrinsic :: iso_c_binding, only: c_int, c_double
      implicit none

      interface
         subroutine VECTOR_PADE(k, n, A, tau,                          &
                   S, gamma, S_star, gamma_star, kappa, num_steps, flag)
            integer,                    intent (in)    :: k
            integer, dimension (:),     intent (in)    :: n
            real,    dimension (:,:),   intent (in)    :: A
            real,                       intent (in)    :: tau
            real,    dimension (:,:,:), intent (out)   :: S, S_star
            real,    dimension (:),     intent (out)   :: gamma,        &
                                                          gamma_star,   &
                                                          kappa
            integer,                    intent (out)   :: num_steps
            integer,                    intent (inout) :: flag
         end subroutine VECTOR_PADE
      end interface

!     bind(C) dummy arguments — all by reference, column-major.
      integer(c_int),                       intent (in)  :: k, nrm
      integer(c_int), dimension(0:k),       intent (in)  :: nvec
      real(c_double), dimension(0:nrm,0:k), intent (in)  :: Amat
      real(c_double),                       intent (in)  :: tau
      real(c_double), dimension(0:nrm),     intent (out) :: Qout
      real(c_double), dimension(0:nrm,1:k), intent (out) :: Pout
      integer(c_int),                       intent (out) :: flag

!     Calgo-766 works in default `real` (single precision in the Sp/
!     build).  We marshal the double-precision caller data down to the
!     single-precision Calgo arrays and back.  The precision floor of
!     this oracle is therefore single precision (~1e-6 to 1e-7); the
!     capture script compares against SharedPade at that tolerance.
      integer                              :: alpha, beta, l, num_steps
      integer, dimension(0:k)              :: nloc
      real,    dimension(0:nrm,0:k)        :: A
      real,    dimension(0:maxval(nvec)+1,0:k,0:k) :: S
      real,    dimension(0:nrm,0:k,0:k)    :: S_star
      real,    dimension(0:k)              :: gamma, gamma_star
      real,    dimension(0:nrm)            :: kappa
      real                                 :: g0

      do beta = 0, k
         nloc(beta) = nvec(beta)
         do l = 0, nrm
            A(l,beta) = real(Amat(l,beta))
         end do
      end do

      flag = 0
      call VECTOR_PADE(k, nloc, A, real(tau),                          &
                S, gamma, S_star, gamma_star, kappa, num_steps, flag)

!     Type-II row-0 result, normalised by gamma_star(0) (Calgo readme:
!     "divide the i'th row of S_star by gamma_star(i)").
      g0 = gamma_star(0)
      Qout = 0.0d0
      Pout = 0.0d0
      do l = 0, nrm
         Qout(l) = real(S_star(l,0,0) / g0, c_double)
      end do
      do beta = 1, k
         do l = 0, nrm
            Pout(l,beta) = real(S_star(l,0,beta) / g0, c_double)
         end do
      end do

      return
      end subroutine vpade_type2
