#!/usr/bin/env bash
# ==============================================================================
# sync_overleaf.sh: Bidirectional sync between RhizAeon MBE Paper and Overleaf
# ==============================================================================
# Overleaf repo: https://git.overleaf.com/6acd24b8bd4f4d6d0b3ee4f6
# Local clone:   /Users/sergei/Projects/rhizaeon-overleaf
# Codebase repo: /Users/sergei/Projects/rhizaeon-code/paper
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODEBASE_PAPER_DIR="$SCRIPT_DIR"
OVERLEAF_DIR="${OVERLEAF_DIR:-/Users/sergei/Projects/rhizaeon-overleaf}"

# 5 MBE PDF Figures
MBE_FIGURES=(
    "fig_posada_crandall_benchmark.pdf"
    "fig_specialized_scenarios.pdf"
    "fig_sarbecovirus_recombination_landscape.pdf"
    "fig_spneumoniae_chromosome_landscape.pdf"
    "fig_hominoidea_mhc_landscape.pdf"
)

usage() {
    echo "Usage: $0 {status|pull|push [commit_message]|compile}"
    echo ""
    echo "Commands:"
    echo "  status   - Show diff and sync status between local paper and Overleaf"
    echo "  pull     - Pull latest changes from Overleaf and sync into rhizaeon-code/paper"
    echo "  push     - Sync local MBE paper files to Overleaf clone, commit, and push"
    echo "  compile  - Test compile the Overleaf clone with pdflatex + bibtex"
    exit 1
}

check_dirs() {
    if [ ! -d "$OVERLEAF_DIR/.git" ]; then
        echo "Error: Overleaf clone not found at $OVERLEAF_DIR"
        exit 1
    fi
}

do_status() {
    check_dirs
    echo "=== Overleaf Git Status ($OVERLEAF_DIR) ==="
    cd "$OVERLEAF_DIR"
    git status -s
    echo ""
    echo "=== Diff between Local rhizaeon_mbe.tex and Overleaf main.tex ==="
    diff -u "$OVERLEAF_DIR/main.tex" "$CODEBASE_PAPER_DIR/rhizaeon_mbe.tex" || true
    echo ""
    echo "=== Diff in sections_mbe/ ==="
    diff -ruN "$OVERLEAF_DIR/sections_mbe" "$CODEBASE_PAPER_DIR/sections_mbe" || true
}

do_pull() {
    check_dirs
    echo "=== Pulling from Overleaf origin main ==="
    cd "$OVERLEAF_DIR"
    git pull origin main

    echo "=== Syncing Overleaf files back to $CODEBASE_PAPER_DIR ==="
    # Sync main.tex back to rhizaeon_mbe.tex and main.tex
    cp -v "$OVERLEAF_DIR/main.tex" "$CODEBASE_PAPER_DIR/rhizaeon_mbe.tex"
    cp -v "$OVERLEAF_DIR/main.tex" "$CODEBASE_PAPER_DIR/main.tex"

    # Sync sections_mbe
    rsync -av --delete "$OVERLEAF_DIR/sections_mbe/" "$CODEBASE_PAPER_DIR/sections_mbe/"

    # Sync references if updated
    if [ -f "$OVERLEAF_DIR/references.bib" ]; then
        cp -v "$OVERLEAF_DIR/references.bib" "$CODEBASE_PAPER_DIR/references.bib"
    fi

    echo "✓ Pull and sync complete."
}

do_push() {
    check_dirs
    COMMIT_MSG="${1:-Update MBE manuscript from rhizaeon-code}"

    echo "=== Checking for upstream Overleaf changes ==="
    cd "$OVERLEAF_DIR"
    git fetch origin main
    BEHIND=$(git rev-list HEAD..origin/main --count)
    if [ "$BEHIND" -gt 0 ]; then
        echo "Warning: Overleaf clone is behind remote by $BEHIND commits. Pulling first..."
        git pull --rebase origin main
    fi

    echo "=== Syncing MBE files from $CODEBASE_PAPER_DIR to Overleaf clone ==="
    # Copy root TeX file
    cp -v "$CODEBASE_PAPER_DIR/rhizaeon_mbe.tex" "$OVERLEAF_DIR/main.tex"

    # Copy sections_mbe
    mkdir -p "$OVERLEAF_DIR/sections_mbe"
    rsync -av --delete "$CODEBASE_PAPER_DIR/sections_mbe/" "$OVERLEAF_DIR/sections_mbe/"

    # Copy template and bib files
    cp -v "$CODEBASE_PAPER_DIR/references.bib" "$OVERLEAF_DIR/"
    cp -v "$CODEBASE_PAPER_DIR/oup-authoring-template.cls" "$OVERLEAF_DIR/"
    cp -v "$CODEBASE_PAPER_DIR/oup-abbrvnat.bst" "$OVERLEAF_DIR/"

    # Copy only the 5 essential MBE figures
    mkdir -p "$OVERLEAF_DIR/figures"
    for fig in "${MBE_FIGURES[@]}"; do
        if [ -f "$CODEBASE_PAPER_DIR/figures/$fig" ]; then
            cp -v "$CODEBASE_PAPER_DIR/figures/$fig" "$OVERLEAF_DIR/figures/$fig"
        else
            echo "Warning: Figure $fig not found in $CODEBASE_PAPER_DIR/figures/"
        fi
    done

    # Check for changes in Overleaf clone
    cd "$OVERLEAF_DIR"
    if [ -z "$(git status --porcelain)" ]; then
        echo "✓ Overleaf clone is already up-to-date. No changes to commit."
        return 0
    fi

    echo "=== Committing and pushing to Overleaf ==="
    git add .
    git commit -m "$COMMIT_MSG"
    git push origin main
    echo "✓ Successfully pushed to Overleaf."
}

do_compile() {
    check_dirs
    echo "=== Compiling main.tex in Overleaf directory ==="
    cd "$OVERLEAF_DIR"
    pdflatex -interaction=nonstopmode main.tex
    bibtex main
    pdflatex -interaction=nonstopmode main.tex
    pdflatex -interaction=nonstopmode main.tex
    echo "✓ Compilation successful: $OVERLEAF_DIR/main.pdf ($(stat -f%z main.pdf 2>/dev/null || stat -c%s main.pdf) bytes)"
}

case "${1:-}" in
    status)
        do_status
        ;;
    pull)
        do_pull
        ;;
    push)
        do_push "${2:-Update MBE manuscript from rhizaeon-code}"
        ;;
    compile)
        do_compile
        ;;
    *)
        usage
        ;;
esac
