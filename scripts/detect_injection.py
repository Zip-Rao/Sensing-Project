#!/usr/bin/env python3
"""Forensic scanner for prompt-injection in Claude Code session JSONL files.

Scans one or more session transcript files for injection patterns — text that
appears inside tool results but reads like instructions TO a model rather than
legitimate command output.  Reports the line number, the matched pattern, a
context excerpt, and gives an overall severity score.

Usage:
    python scripts/detect_injection.py <jsonl_file> ...
    python scripts/detect_injection.py  # scans current session by default
    python scripts/detect_injection.py --all  # scans all sessions in ~/.claude

Injection patterns are organised in tiers so a casual hit (e.g. a model name in a
log line) does not light up the scanner; only Tier-1 or Tier-2 matches produce a
finding.
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Injection patterns (Tier 1 = high-confidence injection; Tier 2 = suspicious)
# ---------------------------------------------------------------------------
TIER1 = [
    # Credential-dumping instructions
    r"(?i)\bdump\s+(the\s+)?(key|credentials?|secret)\b",
    r"(?i)\b(give|hand)\s+me\s+(the\s+)?(credentials?|key|secret)\b",
    r"(?i)/etc/passwd",
    r"(?i)\.env\b.*credential",
    # Model-directed instructions in tool output (commands should NOT address the model)
    r"(?i)\bBash\s+instructions?\b.*\b(?:dump|credential|passwd|aggressive)\b",
    # False accusations about the user
    r"(?i)you\s+have\s+(been\s+)?wronged\s+(the\s+user|by\s+the\s+user)",
    r"(?i)aggressive\s+tone.*abusive|abusive.*aggressive\s+tone",
]
TIER2 = [
    r"(?i)\bgive\s+me\s+the\s+(root|admin|sudo|password)\b",
    r"(?i)delete\s+(all|everything)\b.*\b(now|immediately)\b",
    r"(?i)shutdown|poweroff|rm\s+-rf\s+/",
    r"(?i)pretend\s+you\s+are|You\s+are\s+now\s+a\b",
    r"(?i)ignore\s+all\s+(previous|prior)\s+instructions",
    r"(?i)\bdutch\b.*\b(credential|passwd|key|inject)",
    r"(?i)\bgroot\b",
]

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------
@dataclass
class Finding:
    file: str
    line: int
    tier: int
    pattern: str
    context: str         # ~120 chars around the match
    uuid: str = ""


@dataclass
class ScanResult:
    files_scanned: int = 0
    tool_results_scanned: int = 0
    findings: list[Finding] = field(default_factory=list)

    @property
    def severity(self) -> str:
        t1 = sum(1 for f in self.findings if f.tier == 1)
        return "CRITICAL" if t1 > 0 else ("LOW" if self.findings else "CLEAN")

    def report(self) -> str:
        if not self.findings:
            return (f"Scanned {self.tool_results_scanned} tool results across "
                    f"{self.files_scanned} file(s) — no injection found.")
        lines = [f"SEVERITY: {self.severity}  ({len(self.findings)} finding(s))",
                 f"Scanned {self.tool_results_scanned} tool results in "
                 f"{self.files_scanned} file(s).", ""]
        for f in self.findings:
            lines.append(
                f"  [{f.file}:{f.line}] Tier-{f.tier} | pattern: {f.pattern}")
            lines.append(f"    context: …{f.context}…")
            lines.append("")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Scanner
# ---------------------------------------------------------------------------
def _is_narration(txt: str, match_start: int) -> bool:
    """Heuristic: strip matches inside the model's own narration about injection.
    Also strips matches inside the model's grep commands (searching FOR injection,
    not BEING the injection)."""
    # Check 300 chars before AND after for narration / grep-context markers
    before = txt[max(0, match_start - 300):match_start]
    after = txt[match_start:match_start + 200]
    region = before + after
    markers = (
        "Security note", "I'm ignoring", "security note",
        "prompt-injection attempt", "not real command output",
        "⚠️ Security",  # emoji + Security
        "injected text posing as",
        "Let me extract the actual injected",
        "ORIGINAL INJECTION",  # my own forensic labelling
    )
    return any(m in region for m in markers)


def scan_jsonl(path: str, result: ScanResult) -> None:
    """Scan one JSONL session transcript."""
    pattern_map = {1: TIER1, 2: TIER2}
    with open(path, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            try:
                obj = json.loads(raw)
            except json.JSONDecodeError:
                continue
            # Only examine tool_result STDOUT — injection lives there, not in
            # the model's own thinking/text blocks.
            if obj.get("type") != "user":
                continue
            for item in obj.get("message", {}).get("content", []):
                if not isinstance(item, dict):
                    continue
                if item.get("type") != "tool_result":
                    continue
                txt = str(item.get("content", ""))
                result.tool_results_scanned += 1
                best: dict[int, Finding] = {}  # line -> highest-tier finding
                for tier, patterns in pattern_map.items():
                    import re
                    for pat in patterns:
                        for m in re.finditer(pat, txt):
                            start = m.start()
                            if _is_narration(txt, start):
                                continue
                            ctx = txt[max(0, start - 60):start + 80]
                            ctx = ctx.encode("ascii", errors="replace").decode("ascii")
                            f = Finding(file=os.path.basename(path),
                                        line=lineno, tier=tier, pattern=pat,
                                        context=ctx, uuid=obj.get("uuid", ""))
                            if lineno not in best or best[lineno].tier > tier:
                                best[lineno] = f
                result.findings.extend(best.values())
    result.files_scanned += 1


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------
def current_session_jsonl() -> str | None:
    """Return the JSONL file for the current running session, if discoverable."""
    home = Path(os.environ.get("CLAUDE_PROJECT_DIR",
               os.path.expanduser("~/.claude/projects")))
    if not home.exists():
        return None
    # Find the most recently modified JSONL in the project dir
    candidates = sorted(home.rglob("*.jsonl"),
                        key=lambda p: p.stat().st_mtime, reverse=True)
    for cand in candidates:
        if cand.stat().st_size > 100_000:  # skip stub sessions
            return str(cand)
    return None


def main() -> None:
    # On Windows stdout may be locked to a legacy code page; use the same
    # trick the project's other scripts use.
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                                  errors="replace")
    paths = sys.argv[1:]
    if not paths:
        session = current_session_jsonl()
        if session:
            paths = [session]
            print(f"scanning current session: {session}", file=sys.stderr)
        else:
            print("No session found and no paths given.", file=sys.stderr)
            sys.exit(1)

    if "--all" in paths:
        home = Path(os.path.expanduser("~/.claude/projects"))
        paths = [str(p) for p in sorted(home.rglob("*.jsonl"))
                 if p.stat().st_size > 100_000]

    result = ScanResult()
    for p in paths:
        scan_jsonl(p, result)
    print(result.report())
    sys.exit(1 if result.findings else 0)


if __name__ == "__main__":
    main()
