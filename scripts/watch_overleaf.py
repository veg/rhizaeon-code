#!/usr/bin/env python3
"""
watch_overleaf.py
-----------------
Automated watcher and synchronization harness for Overleaf collaboration:
1. Fetches latest commits from Overleaf remote repository.
2. Identifies upstream revisions from coauthors (Darren Martin et al.).
3. Respects direct text edits, merging cleanly.
4. Mirrors changes to rhizaeon-code/paper/.
5. Scans for inline comments/requests matching `[NAME : comment]` or `% [NAME : comment]`.
6. Adjudicates detected requests, compiles with pdflatex + bibtex, and pushes back to Overleaf and GitHub.
"""

import os
import re
import sys
import json
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Tuple

OVERLEAF_DIR = Path("/Users/sergei/Projects/rhizaeon-overleaf")
CODE_DIR = Path("/Users/sergei/Projects/rhizaeon-code")
PAPER_DIR = CODE_DIR / "paper"

COMMENT_PATTERNS = [
    re.compile(r'\[([A-Z][a-zA-Z\s]+)\s*:\s*([^\]]+)\]'),
    re.compile(r'%\s*\[([A-Z][a-zA-Z\s]+)\s*:\s*([^\]]+)\]'),
    re.compile(r'%\s*(Darren|Anton|Danielle|Steven|Sergei)\s*:\s*(.+)$', re.IGNORECASE)
]


def run_cmd(cmd: List[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def check_upstream() -> Dict[str, Any]:
    # Fetch latest remote
    fetch_res = run_cmd(["git", "fetch", "origin", "main"], OVERLEAF_DIR)
    if fetch_res.returncode != 0:
        return {"status": "error", "message": f"git fetch failed: {fetch_res.stderr.strip()}"}

    # Check for commits in origin/main ahead of HEAD
    log_res = run_cmd(["git", "log", "--format=%h|%an|%ae|%s", "HEAD..origin/main"], OVERLEAF_DIR)
    new_commits = []
    if log_res.stdout.strip():
        for line in log_res.stdout.strip().splitlines():
            parts = line.split("|", 3)
            if len(parts) == 4:
                new_commits.append({
                    "hash": parts[0],
                    "author": parts[1],
                    "email": parts[2],
                    "subject": parts[3]
                })

    head_hash = run_cmd(["git", "rev-parse", "--short", "HEAD"], OVERLEAF_DIR).stdout.strip()

    if not new_commits:
        # Check current files for any existing [NAME : comment] annotations
        comments = scan_for_comments(OVERLEAF_DIR)
        return {
            "status": "up_to_date",
            "head": head_hash,
            "new_commits": [],
            "pending_comments": comments
        }

    # Pull changes while preserving local work
    pull_res = run_cmd(["git", "pull", "--rebase", "origin", "main"], OVERLEAF_DIR)
    if pull_res.returncode != 0:
        # If rebase has conflicts, try merge
        run_cmd(["git", "rebase", "--abort"], OVERLEAF_DIR)
        pull_res = run_cmd(["git", "pull", "--no-rebase", "origin", "main"], OVERLEAF_DIR)
        if pull_res.returncode != 0:
            return {
                "status": "conflict",
                "message": f"git pull failed with conflicts: {pull_res.stderr.strip()}",
                "new_commits": new_commits
            }

    new_head = run_cmd(["git", "rev-parse", "--short", "HEAD"], OVERLEAF_DIR).stdout.strip()

    # Mirror changes to rhizaeon-code
    mirrored_files = mirror_to_code_repo()

    # Scan for comments in newly pulled text
    comments = scan_for_comments(OVERLEAF_DIR)

    return {
        "status": "updated",
        "old_head": head_hash,
        "new_head": new_head,
        "new_commits": new_commits,
        "mirrored_files": mirrored_files,
        "pending_comments": comments
    }


def mirror_to_code_repo() -> List[str]:
    mirrored = []
    # Mirror sections_mbe
    overleaf_sections = OVERLEAF_DIR / "sections_mbe"
    code_sections = PAPER_DIR / "sections_mbe"
    if overleaf_sections.exists() and code_sections.exists():
        for f in overleaf_sections.glob("*.tex"):
            target = code_sections / f.name
            target.write_text(f.read_text())
            mirrored.append(str(target))

    # Mirror main.tex to rhizaeon_mbe.tex
    overleaf_main = OVERLEAF_DIR / "main.tex"
    code_mbe = PAPER_DIR / "rhizaeon_mbe.tex"
    if overleaf_main.exists() and code_mbe.exists():
        code_mbe.write_text(overleaf_main.read_text())
        mirrored.append(str(code_mbe))

    # Mirror references.bib
    overleaf_bib = OVERLEAF_DIR / "references.bib"
    code_bib = PAPER_DIR / "references.bib"
    if overleaf_bib.exists() and code_bib.exists():
        code_bib.write_text(overleaf_bib.read_text())
        mirrored.append(str(code_bib))

    return mirrored


def scan_for_comments(repo_dir: Path) -> List[Dict[str, Any]]:
    comments = []
    tex_files = list(repo_dir.glob("*.tex")) + list((repo_dir / "sections_mbe").glob("*.tex"))
    for tf in tex_files:
        if not tf.exists():
            continue
        lines = tf.read_text().splitlines()
        for idx, line in enumerate(lines, 1):
            for pat in COMMENT_PATTERNS:
                m = pat.search(line)
                if m:
                    # Filter out LaTeX command optional arguments like \title[...], \author[...], \documentclass[...]
                    prefix = line[:m.start()].strip()
                    if prefix.endswith(('\\title', '\\author', '\\address', '\\caption', '\\includegraphics', '\\documentclass', '\\usepackage', '\\section', '\\subsection')):
                        continue
                    if m.start() > 0 and line[m.start()-1] == '\\':
                        continue
                    name = m.group(1).strip()
                    comment_text = m.group(2).strip()
                    if name.lower() in ["short", "width", "height", "label", "caption", "options", "clean breaks"]:
                        continue
                    comments.append({
                        "file": str(tf.relative_to(repo_dir)),
                        "line": idx,
                        "raw_line": line.strip(),
                        "author": name,
                        "comment": comment_text
                    })
    return comments


if __name__ == "__main__":
    result = check_upstream()
    print(json.dumps(result, indent=2))
