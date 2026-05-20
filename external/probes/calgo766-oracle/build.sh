#!/bin/sh
# build.sh — compile ACM Calgo 766 VECTOR_PADE plus the ccall shim into
# shared libraries.  By default this builds the double-precision oracle
# libcalgo766_dp.so, which capture.jl uses; the single-precision build
# libcalgo766.so is kept available for comparison.
#
# Prerequisite: run ./fetch.sh first to populate src-0766/ with the
# Calgo-766 FORTRAN source (gitignored, reproducible from netlib).
#
# Files (all free-form FORTRAN once -ffree-form is forced; the netlib
# .f / .f90 files mix the two layouts, and blas1.f uses free-form `&`
# continuation inside an otherwise fixed-form file, so -ffree-form is
# mandatory and -std=legacy silences the F77-isms):
#   src-0766/Src/Sp/vector_pade.f90  — Calgo 766 core (VECTOR_PADE etc.)
#   src-0766/Src/Sp/linpack.f        — LINPACK routines used by Calgo
#   src-0766/Src/Sp/blas1.f          — BLAS-1 routines used by LINPACK
#   wrapper.f90                      — working_area module + bind(C) shim
#
# PRECISION
# ---------
# The netlib TOMS mirror (https://www.netlib.org/toms/766) ships a single
# `Sp/` source tree only — there is no `Dp/` double-precision variant.
# But the Calgo 766 FORTRAN is fully precision-agnostic: every floating
# declaration is a bare `real` (no `kind` parameter, no `selected_real_kind`,
# no `real*8`), and the lone numeric literal in linpack.f (`-1.0e0`) is a
# default-`real`-kind literal.  Compiling with
#   -fdefault-real-8 -fdefault-double-8
# therefore promotes the entire algorithm to IEEE double precision with no
# source edits.  The `bind(C)` ABI in wrapper.f90 is expressed via
# `iso_c_binding` `c_double` / `c_int` named constants, which are NOT
# affected by -fdefault-real-8 (they always denote the C `double` / `int`
# kinds), so the ccall ABI in capture.jl stays Float64 / int32 unchanged.
# The `real(...)` marshalling conversions in wrapper.f90 become double->
# double no-ops in the Dp build.
#
# Accuracy floor:
#   libcalgo766.so      single precision — oracle floor ~1e-6..1e-7
#   libcalgo766_dp.so   double precision — oracle floor ~1e-15 (machine eps)
set -e
cd "$(dirname "$0")"

SRC=src-0766/Src/Sp
if [ ! -f "$SRC/vector_pade.f90" ]; then
  echo "build.sh: $SRC/vector_pade.f90 missing — run ./fetch.sh first." >&2
  exit 1
fi

BASE_FFLAGS="-ffree-form -std=legacy -O2 -fPIC"

# wrapper.f90 defines the module `working_area_VECTOR_PADE` that
# vector_pade.f90 USEs.  The precision of that module's variables is baked
# into the generated working_area_vector_pade.mod, so the single- and
# double-precision builds MUST NOT share a .mod file (or a stale .mod in
# the current directory).  Each build therefore compiles in its own private
# directory with -J pointing there; nothing is written to the probe dir.
build_variant() {
    # $1 = output .so name   $2 = private build dir   $3+ = extra FFLAGS
    out_so="$1"; bdir="$2"; shift 2
    rm -rf "$bdir"; mkdir -p "$bdir"
    fflags="$BASE_FFLAGS $* -J$bdir"
    # wrapper.f90 first: it produces working_area_vector_pade.mod into $bdir.
    gfortran $fflags -c wrapper.f90            -o "$bdir/wrapper.o"
    gfortran $fflags -c "$SRC/vector_pade.f90" -o "$bdir/vector_pade.o"
    gfortran $fflags -c "$SRC/linpack.f"       -o "$bdir/linpack.o"
    gfortran $fflags -c "$SRC/blas1.f"         -o "$bdir/blas1.o"
    gfortran -shared -o "$out_so" \
        "$bdir/wrapper.o" "$bdir/vector_pade.o" \
        "$bdir/linpack.o" "$bdir/blas1.o"
    rm -rf "$bdir"
}

# --- double-precision build (the oracle capture.jl uses) ---------------------
# -fdefault-real-8 promotes bare `real` to 8 bytes; -fdefault-double-8 keeps
# `double precision` consistent.  c_double / c_int are unaffected.
build_variant libcalgo766_dp.so .dp-build -fdefault-real-8 -fdefault-double-8
echo "build.sh: wrote $(pwd)/libcalgo766_dp.so  (double precision)"

# --- single-precision build (kept for comparison) ---------------------------
build_variant libcalgo766.so .sp-build
echo "build.sh: wrote $(pwd)/libcalgo766.so     (single precision)"
