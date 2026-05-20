#!/bin/sh
# build.sh — compile ACM Calgo 766 VECTOR_PADE plus the ccall shim into
# a shared library libcalgo766.so.
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
# Calgo 766 in the Sp/ build uses default `real` = single precision, so
# the oracle's accuracy floor is ~1e-6..1e-7 (see wrapper.f90 header).
set -e
cd "$(dirname "$0")"

SRC=src-0766/Src/Sp
if [ ! -f "$SRC/vector_pade.f90" ]; then
  echo "build.sh: $SRC/vector_pade.f90 missing — run ./fetch.sh first." >&2
  exit 1
fi

FFLAGS="-ffree-form -std=legacy -O2 -fPIC"

# wrapper.f90 must be compiled first: it defines working_area_VECTOR_PADE,
# the module vector_pade.f90 USEs.
gfortran $FFLAGS -c wrapper.f90              -o wrapper.o
gfortran $FFLAGS -c "$SRC/vector_pade.f90"   -o vector_pade.o
gfortran $FFLAGS -c "$SRC/linpack.f"         -o linpack.o
gfortran $FFLAGS -c "$SRC/blas1.f"           -o blas1.o

gfortran -shared -o libcalgo766.so \
    wrapper.o vector_pade.o linpack.o blas1.o

echo "build.sh: wrote $(pwd)/libcalgo766.so"
