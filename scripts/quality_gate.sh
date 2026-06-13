#!/usr/bin/env bash
# =============================================================================
# scripts/quality_gate.sh — the single, STRICTLY-SERIAL local quality-gate runner.
#
# Bead padetaylor-krgy.10 (epic padetaylor-krgy, EFFICIENCY/§Runner of
# docs/test_corpus/03_hardening_methodology.md:171-174). The CAPSTONE of the
# test-hardening sweep: it does not add a gate, it WIRES the gates the sweep
# built into one command a human can run before committing or pushing.
#
# WHY THIS EXISTS. There is NO GitHub CI here (CLAUDE.md Rule 11) — quality
# gates only protect the repo if they are trivially runnable locally, and only
# help a human who has no CI log to read if their PASS/FAIL is interpreted FOR
# them. This script is that one entry point: it sequences the gates, parses each
# one's result, and prints a verdict that reads as GREEN when the suite is
# healthy — crucially WITHOUT alarming the human over the two facts a naive
# reader misreads as failure (see "EXPECTED-NOISE" below).
#
# THE SERIAL INVARIANT (CLAUDE.md Rule 7 — load-bearing). Parallel `julia`
# corrupts the precompile cache, so this repo runs ONE Julia process at a time.
# This runner honours that BY CONSTRUCTION: it is plain bash that launches each
# gate as a SEPARATE child `julia` invocation, in sequence, each fully exiting
# before the next is launched. There is no `&`, no job control, no background
# subshell anywhere below — every gate is a foreground call whose exit code we
# capture before moving on. ReTestItems.jl was CONSIDERED AND REJECTED for the
# suite precisely because its worker-process parallelism violates Rule 7
# (methodology doc; bead NOTES, confirmed 3-0). bash sequencing separate
# processes is the inherently Rule-7-safe shape.
#
# EXPECTED-NOISE (the two facts a passing run must NOT look alarming over):
#   1. The full Pkg.test() suite is GREEN at N pass / 3 BROKEN / 0 fail. The
#      3 broken are INTENTIONAL @test_broken markers: 1 auto-flip marker for a
#      known open bug (v1ub×1, CFail.1d) + 2 deferred-feature markers
#      (cm-n2 collision in corpus_vector_polefield_test, pi2-tritronquee in
#      corpus_painleve_rational_test). [Was 10; on 2026-06-13 padetaylor-53tu's
#      3 markers (CBvx.4.2, CBvx.4.4, CBV.7) AND padetaylor-q0yq's 1 marker
#      (CBV.9) were FIXED and flipped to @test_throws; then padetaylor-61um's
#      3 markers (CWD.5, CBr.3, CPN.7.3) were REFRAMED 6→3 (winding_delta is
#      correct; the fix is BranchTracker's grazing guard — ADR-0031).  NB the
#      conditional kkg_pi2_figure_test markers contribute 0 when the Stage-B
#      march succeeds, as in the canonical run.] A nonzero broken count is EXPECTED, not a regression — each
#      auto-flip marker flips to "Unexpected Pass" the day its bug is fixed
#      (that IS the signal to remove the marker). INVESTIGATE ONLY FAILs.
#      (bd memory corpus-v2-expected-broken-count.) So this runner gates on the
#      `N failed` count it parses out of the summary, NOT on "broken == 0".
#   2. The mutation gate carries 2 TRACKED survivors under bead padetaylor-98pe
#      (the 2nd-order up-resync index slip + the SVD rank-threshold mutant — both
#      catalogued :survived). A tracked survivor is NOT a failure; the gate exits
#      0 when every mutant meets its catalogued expectation. (test/mutation/README.md.)
# The runner SURFACES both facts in its summary so the human is not alarmed.
#
# TIERS (argument-selected; default `full`):
#   fast  — static-analysis gate only (test/quality_test.jl: Aqua + JET +
#           ExplicitImports), ~1.5 min cold. Pre-commit. NOTE: this is the
#           STATIC-ANALYSIS subset only. The type-stability gate is deliberately
#           NOT in `fast`: it depends on AllocCheck, which is a test-target dep
#           NOT resolvable in the plain `--project=.` env (it lives only inside
#           the Pkg.test sandbox), and reconstructing that sandbox by hand would
#           be exactly the brittle shortcut a senior does not ship (Rule 9). Run
#           `full` to exercise type-stability. quality_test.jl, by contrast, is
#           self-contained and its deps stack cleanly outside the sandbox —
#           verified standalone GREEN, krgy.10.
#   full  — the authoritative correctness suite: `Pkg.test()`. ~24 min. The
#           always-on gate; includes ALL in-suite gates (static analysis,
#           property, metamorphic, certified-oracle, differential, convergence,
#           accuracy-ledger, oracle-provenance, type-stability — see
#           test/runtests.jl). The pre-push gate.
#   deep  — `full` THEN the three PERIODIC out-of-suite tools in sequence:
#           mutation gate → coverage snapshot → perf gate. End-of-phase only.
#           ~24 + ~5 (mutation) + ~24 (coverage regen) + ~2 (perf) ≈ 55 min, and
#           the perf gate is meaningful ONLY on an IDLE box with threads pinned
#           (it reports its thread config and refuses to regen otherwise —
#           benchmark/README.md idle-box discipline).
#
# Each gate is its own Julia child; full `deep` runs FOUR julia processes back
# to back, never two at once.
#
# USAGE:
#   scripts/quality_gate.sh [fast|full|deep] [--dry-run]
#   scripts/quality_gate.sh --help | --list
# `--dry-run` prints the exact, ordered commands a tier would run (and exits 0)
# without launching any julia — use it to inspect/verify the plumbing.
#
# NOT a Pkg.test include (it IS the orchestrator); adds no dependency; touches
# no src/. See CLAUDE.md "Practical guidance" for the pointer.
# =============================================================================
set -u -o pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
JULIA="${JULIA:-julia}"

usage() {
  sed -n '2,/^# ===.*===$/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

# ── Per-tier command plans (printed verbatim by --dry-run; run in order) ───────
# Each element is a "<label>|<command>" the runner executes as ONE julia child.
fast_plan() {
  echo "static-analysis (quality_test.jl)|$JULIA --project=$REPO $REPO/test/quality_test.jl"
}
full_plan() {
  echo "correctness suite (Pkg.test)|$JULIA --project=$REPO -e 'using Pkg; Pkg.test()'"
}
deep_plan() {
  full_plan
  echo "mutation gate|$JULIA --project=$REPO $REPO/test/mutation/run_mutation_gate.jl"
  echo "coverage snapshot|$JULIA --project=$REPO $REPO/test/coverage/run_coverage.jl --run"
  echo "perf gate (IDLE box only)|$JULIA --project=$REPO/benchmark $REPO/benchmark/run_perf_gate.jl"
}

plan_for() {
  case "$1" in
    fast) fast_plan ;;
    full) full_plan ;;
    deep) deep_plan ;;
    *) echo "quality_gate: unknown tier '$1' (want fast|full|deep)" >&2; return 2 ;;
  esac
}

# ── Test-summary table parse (POSITIONAL — load-bearing) ──────────────────────
# Julia's top-level `Test Summary:` prints a two-row table: a HEADER row naming
# the columns (`Pass`/`Fail`/`Broken`/`Total`/`Time`, in that order, with `Fail`
# and `Broken` present ONLY when nonzero) and a DATA row whose integers align
# UNDER those headers by COLUMN INDEX, not adjacency:
#     Test Summary: | Pass  Broken  Total      Time
#     PadeTaylor.jl | 8859      10   8869  18m30.0s
# So we must map the i-th numeric token in the data row to the i-th column word
# in the header row — a word-adjacency regex (`[0-9]+ *Broken`) FAILS here
# because the number and its header word are on different lines. Emits the three
# counts as "pass fail broken" (missing ⇒ 0); fail="?" only if no table is found.
_suite_counts() {
  awk '
    /Test Summary:/ && /Pass/ { split($0, h, "|"); n=split(h[2], col, /[ \t]+/);
      idx=0; for (i=1;i<=n;i++){ if(col[i]!=""){ idx++; name[idx]=col[i] } }
      ncol=idx; getline; split($0, d, "|"); m=split(d[2], num, /[ \t]+/);
      j=0; for (i=1;i<=m;i++){ if(num[i]!=""){ j++; val[name[j]]=num[i] } }
      p=(val["Pass"]==""?0:val["Pass"]); f=(val["Fail"]==""?"?":val["Fail"]);
      b=(val["Broken"]==""?0:val["Broken"]); print p, f, b; found=1; exit }
    END { if (!found) print "?", "?", "0" }
  ' "$1"
}

# ── Result interpretation ──────────────────────────────────────────────────────
# A captured log + exit code from one gate → a PASS/FAIL verdict line. The full
# suite needs special handling: its exit code already reflects FAILs, but we ALSO
# parse the count table (above) so we can echo the EXPECTED-NOISE reassurance
# instead of letting a human stare at "10 broken" and panic.
interpret() {
  local label="$1" rc="$2" log="$3"
  if [[ "$label" == correctness* ]]; then
    local passed failed broken
    read -r passed failed broken < <(_suite_counts "$log")
    # No `Fail` column in the table ⇒ zero failures (Julia omits it when zero).
    [[ "$failed" == "?" ]] && grep -q "Test Summary:" "$log" && failed=0
    if [[ "$rc" -eq 0 && "$failed" == "0" ]]; then
      echo "  PASS  correctness suite — $passed pass / $broken broken / 0 fail"
      echo "        ($broken broken is EXPECTED: intentional @test_broken markers —"
      echo "         v1ub auto-flip + 2 deferred-feature (cm-n2, pi2-tritronquee)."
      echo "         q0yq/53tu/61um were FIXED. NOT a regression; investigate only"
      echo "         FAILs. bd memory corpus-v2-expected-broken-count.)"
      return 0
    fi
    echo "  FAIL  correctness suite — $failed failed (exit $rc). INVESTIGATE: real"
    echo "        FAILs (not the $broken broken markers). Re-read the log above."
    return 1
  fi
  if [[ "$label" == mutation* ]]; then
    [[ "$rc" -eq 0 ]] && { echo "  PASS  mutation gate — exit 0 (every mutant met its catalogued"; \
      echo "        expectation; the 2 tracked survivors under bead padetaylor-98pe are"; \
      echo "        EXPECTED, not failures — test/mutation/README.md)."; return 0; }
    echo "  FAIL  mutation gate — exit $rc: an UNEXPECTED outcome (a once-killed mutant"
    echo "        now survives ⇒ weakened test, OR a tracked survivor got killed ⇒"
    echo "        update the catalogue). See the SURVIVORS block above."
    return 1
  fi
  if [[ "$label" == coverage* ]]; then
    echo "  $([[ "$rc" -eq 0 ]] && echo INFO || echo FAIL)  coverage snapshot — exit $rc"
    echo "        (a MEASUREMENT, not a pass/fail gate; read the per-file % above and"
    echo "         triage cold lines per docs/test_corpus/coverage_snapshot.md)."
    return "$rc"
  fi
  # generic (perf + fast static-analysis)
  if [[ "$rc" -eq 0 ]]; then echo "  PASS  $label — exit 0."; return 0; fi
  echo "  FAIL  $label — exit $rc. Re-read the log above."
  return 1
}

# ── Driver: launch each gate SERIALLY, capture, interpret, accumulate ─────────
run_tier() {
  local tier="$1" dry="$2" plan rc=0 step=0
  plan="$(plan_for "$tier")" || return $?
  local n; n=$(printf '%s\n' "$plan" | grep -c .)
  echo "=== quality_gate: tier '$tier' — $n gate(s), STRICTLY SERIAL (Rule 7) ==="
  if [[ "$dry" == "1" ]]; then
    echo "--- DRY RUN: the ordered commands this tier would run (no julia launched) ---"
    while IFS='|' read -r label cmd; do step=$((step+1)); printf '  [%d/%d] %-32s %s\n' "$step" "$n" "$label" "$cmd"; done <<< "$plan"
    return 0
  fi
  local summary=() overall=0 log
  while IFS='|' read -r label cmd; do
    step=$((step+1))
    echo; echo ">>> [$step/$n] $label"; echo ">>> $cmd"
    log="$(mktemp)"
    # ONE foreground julia child; fully exits before the loop advances (Rule 7).
    ( cd "$REPO" && eval "$cmd" ) 2>&1 | tee "$log"; rc=${PIPESTATUS[0]}
    local verdict; verdict="$(interpret "$label" "$rc" "$log")" || overall=1
    echo "$verdict"; summary+=("$verdict")
    rm -f "$log"
  done <<< "$plan"
  echo; echo "================ quality_gate SUMMARY (tier '$tier') ================"
  printf '%s\n' "${summary[@]}"
  echo "===================================================================="
  [[ "$overall" -eq 0 ]] && echo "VERDICT: GREEN — all gates passed (expected broken/survivors noted above)." \
                         || echo "VERDICT: RED — a gate FAILED; see the FAIL line(s) above."
  return "$overall"
}

# ── Argument parsing ───────────────────────────────────────────────────────────
TIER="full"; DRY="0"
for a in "$@"; do
  case "$a" in
    -h|--help) usage; exit 0 ;;
    --list) echo "tiers: fast | full (default) | deep"; exit 0 ;;
    --dry-run) DRY="1" ;;
    fast|full|deep) TIER="$a" ;;
    *) echo "quality_gate: unknown arg '$a' (see --help)" >&2; exit 2 ;;
  esac
done
run_tier "$TIER" "$DRY"
