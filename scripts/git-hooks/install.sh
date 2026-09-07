#!/usr/bin/env bash
#
# install.sh — point this clone's git hooks at the repo-tracked ones.
#
# `core.hooksPath` is per-clone configuration and is NOT carried by a clone
# or a pull, so every machine that works on this repo runs this once.  The
# hooks themselves live in the repo, so they stay identical across machines.
#
#   ./scripts/git-hooks/install.sh
#
# Undo with:  git config --unset core.hooksPath

set -euo pipefail
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"
chmod +x scripts/git-hooks/pre-commit
git config core.hooksPath scripts/git-hooks
echo "installed: core.hooksPath -> scripts/git-hooks"
echo "hooks active:"
for h in scripts/git-hooks/*; do
    case "$h" in */install.sh) continue ;; esac
    echo "  - $(basename "$h")"
done
