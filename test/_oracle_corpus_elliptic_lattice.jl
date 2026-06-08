# test/_oracle_corpus_elliptic_lattice.jl
#
# Pinned mpmath oracle constants for the ELLIPTIC-LATTICE corpus (bead
# padetaylor-dctk, epic padetaylor-25og; plan Family C,
# docs/test_corpus/02_corpus_extension_plan.md:204-228).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/elliptic-lattices/capture.py  (sympy 1.14
#   + mpmath 1.3, mp.dps = 50).  Every constant below is a verbatim copy of
#   that script's labelled output.  Before any value is emitted the script
#   GATES on the recast-ODE residual AND a closed-form invariant:
#     * CEl.1: u''=6u^2-1/2 residual ≤1e-38 and (P')²=4P³-g₂P-g₃ ≤1e-38;
#     * CEl.2: cn''/dn'' recast residuals ≤1e-38 and the first-order
#       identities (cn')²=(1-cn²)(1-m+m cn²), (dn')²=(1-dn²)(dn²-1+m);
#     * CEl.3: the d=3 SYSTEM residuals sn'=cn dn / cn'=-sn dn / dn'=-m sn cn
#       ≤1e-38, the invariants sn²+cn²=1, dn²+m sn²=1 ≤1e-40, AND the
#       far-side algebraic identity (sn,cn,dn)=(-i/√k,-√((1+k)/k),-√(1+k))
#       at z=i·3K'/2 ≤1e-40.
#   To update: re-run `python3 capture.py` and re-pin.
#
# Ground-truth recipe: plan Family C (md:204-228); DLMF 22.x (Jacobi
# cn/dn/sn) + 23.x (Weierstrass ℘).  The lemniscatic lattice g₂=1, g₃=0 is
# the REAL-root case (e = {1/2, 0, −1/2}) where the validated ℘-via-Jacobi-sn
# recipe is dps=50-clean (plan discovery #5).  The Jacobi cases reuse mpmath
# ellipfun directly.
#
# Exact ALGEBRAIC constants (the −i/√k far-side targets, the ±400 bare pole)
# are written INLINE / computed in-test, not here — only transcendental pins
# live in this file.

# ===========================================================================
# CEl.1  lemniscatic ℘ (SQUARE lattice):  u'' = 6u² − 1/2, u = ℘(z−1; 1, 0).
# IC at z₀=0; double pole shifted to z=1 (parallel to the v1 CPB.1 keystone).
# ===========================================================================
const CEL1_U0       = 1.05083979104023701837944993034901171657327564
const CEL1_UP0      = 1.89493523185581089049547932437906567542905606
const CEL1_U_0p5    = 4.01251302709622740374047352926338115141493875
const CEL1_U_1p05   = 400.000125000013020833959334964665893096310149   # PAST pole
const CEL1_UP_1p05  = -15999.9949999984374998747996714320113755998357
const CEL1_U_1p4    = 6.25800341400562398250212279113182779665756603    # far past
# The LATTICE-DISCRIMINATING signal (plan md:210): at z=1.05 the bare double
# pole 1/(z−1)² = 400 exactly; the SQUARE-lattice correction above 400 is
# ~1.25e-4 (lemniscatic) vs ~4.46e-7 for the equianharmonic 60° rhombus — an
# ~280× contrast at FIXED (double) pole order.  Both pinned so the test
# discriminates the lattice GEOMETRY, not just the local pole.
const CEL1_LATTICE_CORR_1p05 = 0.000125000013020833959334964665893096310149088 # square
const CEL1_EQUI_CORR_1p05    = 0.000000446428571466898057302094343353962712214 # rhombus

# ===========================================================================
# CEl.2  Jacobi cn, dn (scalar 2nd-order; complex bridge through iK').
# cn: u''=(2m−1)u−2m u³, IC (1,0);  dn: u''=(2−m)u−2u³, IC (1,0).  m=0.36
# (k=0.6).  Simple pole at z=iK', residues −i/k (cn), −i (dn).  cn and dn are
# REAL on the imaginary axis between lattice points; they flip sign across the
# pole (PRE at z=iK'/2 positive, POST at z=i·3K'/2 negative).
# ===========================================================================
const CEL2_K  = 1.75075380291575252897522604601214825576745916   # K(0.36)
const CEL2_KP = 1.99530277766472938768621133937243734938196807   # K'(0.36)=K(0.64)
# cn / dn at z = iK'/2 (BEFORE the pole) and z = i·3K'/2 (PAST the pole).
const CEL2_CN_PRE  = 1.63299316185545208371198041146740065790424345   # cn(iK'/2)
const CEL2_DN_PRE  = 1.26491106406735174128005447004088027993228895   # dn(iK'/2)
const CEL2_CN_POST = -1.63299316185545219202214369954929972759330398  # cn(i·3K'/2)
const CEL2_DN_POST = -1.26491106406735179161806950681053942268086180  # dn(i·3K'/2)
# Residue imaginary parts (k=0.6): cn ~ (−i/k)/(z−iK'), dn ~ (−i)/(z−iK').
const CEL2_RES_CN_IM = -1.66666666666666666666666666666666666666666667  # = −1/k
const CEL2_RES_DN_IM = -1.0                                            # = −1

# ===========================================================================
# CEl.3  Jacobi TRIPLE (first-order d=3 vector — the shared-Q oracle).
# sn'=cn dn, cn'=−sn dn, dn'=−m sn cn, IC (0,1,1).  m=1/2 (self-dual K=K',
# square lattice, links to CEl.1).  Three components share ONE identical pole
# lattice (simple poles at iK').  Far-side EXACT algebraic at z=i·3K'/2 (the
# half-period past iK'): (sn,cn,dn)=(−i/√k, −√((1+k)/k), −√(1+k)), k=√m.
# ===========================================================================
const CEL3_K  = 1.85407467730137191843385034719526004621759882   # K(1/2)=K'(1/2)
const CEL3_KP = 1.85407467730137191843385034719526004621759882   # self-dual
# PRE-pole values at z = iK'/2 (sn imaginary, cn/dn real).
const CEL3_SN_PRE_IM = 1.18920711500272102387408713795939172605272856   # sn(iK'/2)
const CEL3_CN_PRE    = 1.55377397403003727455323039903451480166733731   # cn(iK'/2)
const CEL3_DN_PRE    = 1.30656296487637650835904040142881843051511725   # dn(iK'/2)
# POST-pole values at z = i·3K'/2 (the mpmath values, with the Float64 K'
# cast the test's zspan also uses — so a ~1e-14 agreement is the realistic
# bar against this pin; the EXACT algebraic targets below are the cleaner
# cross-check the test also pins).
const CEL3_SN_POST_IM = -1.18920711500272142063450663471457724180195010 # sn(i·3K'/2)
const CEL3_CN_POST    = -1.55377397403003757822050871829684773672574408 # cn(i·3K'/2)
const CEL3_DN_POST    = -1.30656296487637668892068438681797588338105606 # dn(i·3K'/2)
# Far-side EXACT algebraic (k=√(1/2)=2^(−1/2)); computed in-test, pinned here
# only as the gate's reference (capture.py GATE C confirms they equal the
# full-precision sn/cn/dn at i·3K'/2 to ≤1e-40).
const CEL3_SN_ALG_IM = -1.18920711500272106671749997056047591529297209 # −1/√k
const CEL3_CN_ALG    = -1.55377397403003730734415895306314694816458350 # −√((1+k)/k)
const CEL3_DN_ALG    = -1.30656296487637652785664317342718715358376119 # −√(1+k)
