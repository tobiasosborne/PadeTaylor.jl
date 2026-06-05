#!/usr/bin/env python3
"""Deterministically render the Phase-2 corpus catalogue from the workflow JSON.

The ground-truth RECIPES are load-bearing (they regenerate the reference values a
test asserts against), so they must be reproduced byte-exactly — this script copies
fields verbatim rather than risk an LLM paraphrase.  Re-run to regenerate:

    python3 docs/test_corpus/_generate_catalogue.py <workflow.output.json>

Writes docs/test_corpus/01_corpus_catalogue.md and docs/test_corpus/candidates.json.
"""
import json, sys, os

src = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(src))["result"]
territories = d["territories"]
con = d["consolidated"]

# persist the raw candidate records for provenance / subagent reference
json.dump(d, open(os.path.join(here, "candidates.json"), "w"), indent=1)

out = []
W = out.append
W("# PadeTaylor.jl Ground-Truth Corpus Catalogue\n")
W("Phase-2 deliverable of epic `padetaylor-p1v0` (comprehensive ground-truth test "
  "corpus). Gathered by 9 read-only Sonnet scout territories + a consolidation pass, "
  "under a strict anti-hallucination contract: every entry's ground truth is a "
  "**closed-form** solution, a **regenerable recipe** in locally-available tooling "
  "(Mathematica `wolframscript`, Python `mpmath`/`scipy`/`sympy`, Julia "
  "`DifferentialEquations`/`HypergeometricFunctions`/`Nemo`), or a **precisely-cited "
  "published value** — never a fabricated number. Recipes are reproduced verbatim "
  "from `candidates.json` by `_generate_catalogue.py`.\n")
W("> Companion: `00_capability_map.md` (the 22 capability buckets B1–B22 + 15 "
  "highest-value gaps this corpus targets).\n")

# ---- coverage matrix -------------------------------------------------------
W("## Bucket coverage after the sweep\n")
W("| Bucket | Coverage | Best candidates | Remaining hole |")
W("|---|---|---|---|")
order = {"strong-regenerable": 0, "partial": 1, "strong-citation-only": 2, "still-uncovered": 3}
for b in sorted(con["bucket_coverage"], key=lambda x: x["bucket"]):
    best = ", ".join(f"`{c}`" for c in b.get("best_candidates", []))
    hole = (b.get("remaining_hole") or "—").replace("\n", " ").replace("|", "\\|")
    W(f"| {b['bucket']} | {b['now_covered']} | {best} | {hole} |")
W("")

W("## Still uncovered — source manually\n")
for g in con["still_uncovered"]:
    W(f"- {g}")
W("")

W("## Top-priority oracles (implement first)\n")
for i, g in enumerate(con["top_priority_oracles"], 1):
    W(f"{i}. `{g}`")
W("")

# ---- consolidated index ----------------------------------------------------
W("## Consolidated corpus index\n")
W("`R` = regenerable (closed-form or local recipe); `C` = citation-only.\n")
W("| Pri | R/C | Buckets | id | gt_kind | one-line |")
W("|---|---|---|---|---|---|")
for e in sorted(con["corpus_index"], key=lambda x: (x["priority"], x["buckets"])):
    rc = "R" if e.get("regenerable") else "C"
    ol = e["one_line"].replace("\n", " ").replace("|", "\\|")
    W(f"| {e['priority']} | {rc} | {','.join(e['buckets'])} | `{e['id']}` | {e['gt_kind']} | {ol} |")
W("")

# ---- full per-candidate detail --------------------------------------------
W("## Full candidate records\n")
W("Every field verbatim from the scout output. `how_to_compute` is the regeneration "
  "recipe / source locator a test's oracle must use.\n")
FIELDS = ["name", "ode", "ic_bc", "domain", "closed_form", "pole_structure",
          "gt_kind", "how_to_compute", "precision", "targets_buckets", "regime",
          "citation", "confidence", "verify_note"]
total = 0
for t in territories:
    W(f"### Territory: {t['territory'][:80]}\n")
    if t.get("territory_notes"):
        W(f"*Scout notes:* {t['territory_notes']}\n")
    for c in t.get("candidates", []):
        total += 1
        W(f"#### `{c.get('id','?')}`\n")
        for f in FIELDS:
            if f in c and c[f] not in (None, "", []):
                v = c[f]
                if isinstance(v, list):
                    v = ", ".join(map(str, v))
                v = str(v).replace("\n", " ")
                W(f"- **{f}**: {v}")
        W("")

W(f"\n---\n*{total} candidate records across {len(territories)} territories; "
  f"{len(con['corpus_index'])} consolidated entries.*")

open(os.path.join(here, "01_corpus_catalogue.md"), "w").write("\n".join(out) + "\n")
print(f"wrote 01_corpus_catalogue.md ({len(out)} lines), candidates.json, "
      f"{total} candidate records")
