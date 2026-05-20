#!/bin/sh
# fetch.sh — reproducibly acquire the ACM Calgo 766 FORTRAN source.
#
# Calgo 766: Cabay, Jones & Labahn, "Algorithm 766: Experiments with a
# Weakly Stable Algorithm for Computing Pade-Hermite and Simultaneous
# Pade Approximants", ACM TOMS 23(1), 91-110, March 1997,
# DOI 10.1145/244768.244790.
#
# The journal paper is paywalled; the FORTRAN source is openly
# redistributable.  The ACM Calgo archive (https://calgo.acm.org/0766.zip)
# now 404s, but the netlib TOMS mirror still serves it as a gzipped Unix
# shell archive at https://www.netlib.org/toms/766.
#
# This script downloads that archive, gunzips it, and runs the embedded
# shar to populate src-0766/ with:
#   Doc/readme            Drivers/Sp/{driver1,driver2}.f90  Drivers/Sp/data2
#   Src/Sp/vector_pade.f90  Src/Sp/linpack.f  Src/Sp/blas1.f  ...
#
# src-0766/ is gitignored (see .gitignore); this script + a network
# connection reproduces it.
set -e
cd "$(dirname "$0")"

URL="https://www.netlib.org/toms/766"

echo "fetch.sh: downloading Calgo 766 from $URL"
curl -fsSL -o toms766.gz "$URL"

echo "fetch.sh: gunzipping (netlib serves it gzip-compressed)"
gunzip -c toms766.gz > toms766.shar

rm -rf src-0766
mkdir -p src-0766
cp toms766.shar src-0766/766.shar
# The shar has 4 banner lines beginning 'C ' before the '#! /bin/sh'
# line; `sh` ignores them with harmless 'C: not found' warnings.
( cd src-0766 && sh 766.shar >/dev/null 2>&1 || true )

if [ ! -f src-0766/Src/Sp/vector_pade.f90 ]; then
  echo "fetch.sh: ERROR — extraction failed, vector_pade.f90 not found." >&2
  exit 1
fi

rm -f toms766.gz toms766.shar
echo "fetch.sh: Calgo 766 source extracted to src-0766/"
