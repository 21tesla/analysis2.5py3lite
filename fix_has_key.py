#!/usr/bin/env python
"""Migrate dict.has_key(x) → x in dict across the codebase.

Handles:
  - obj.has_key(key)      → key in obj
  - not obj.has_key(key)  → key not in obj
  - Complex left-hand expressions: self.x, d['key'], module.attr, etc.

Safety: only rewrites lines where has_key appears as a method call
with a simple (non-nested-paren) argument.
"""
import re
import sys
from pathlib import Path

# Match: <expr>.has_key(<arg>)
#   <expr>: word chars, dots, attribute access, subscript [..]
#   <arg>:  anything that doesn't contain a ) 
# We capture the expr and arg separately.
HAS_KEY_RE = re.compile(
    r"([a-zA-Z_]\w*(?:\.[a-zA-Z_]\w*|\[[^\]]+\]|\"[^\"]*\"|'[^']*')*)"
    r"(?:\.has_key\("
    r"([^\)]+)"
    r"\))"
)

def fix_line(line: str) -> str:
    """Fix a single line, rewriting has_key calls to `in` expressions."""
    match = HAS_KEY_RE.search(line)
    if not match:
        return line

    obj = match.group(1)
    arg = match.group(2)
    start, end = match.span()

    # Check if preceded by "not " (with optional whitespace)
    prefix = line[:start].rstrip()
    if prefix.endswith("not"):
        replacement = f"{arg} not in {obj}"
        # Remove the "not " prefix
        prefix = prefix[:-3].rstrip()
        if prefix.endswith(" ") and not prefix.endswith("  "):
            pass  # keep one space
        result = prefix + " " + replacement + line[end:]
    else:
        replacement = f"{arg} in {obj}"
        result = line[:start] + replacement + line[end:]

    # Recurse for multiple has_key on the same line
    if ".has_key(" in result:
        return fix_line(result)
    return result


def process_file(path: Path) -> bool:
    """Process one file, return True if modified."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if ".has_key(" not in text:
        return False

    new_lines = []
    for line in text.splitlines(keepends=True):
        stripped = line.rstrip("\n")
        newline = "\n" if line.endswith("\n") else ""
        if ".has_key(" in stripped:
            stripped = fix_line(stripped)
        new_lines.append(stripped + newline)

    new_text = "".join(new_lines)
    if new_text != text:
        path.write_text(new_text, encoding="utf-8")
        return True
    return False


def main(root: str):
    root_path = Path(root)
    changed = 0
    total = 0
    for py_file in sorted(root_path.rglob("*.py")):
        if "__pycache__" in str(py_file):
            continue
        total += 1
        if process_file(py_file):
            changed += 1
            print(f"  fixed: {py_file}")
    print(f"\n{changed}/{total} files modified.")

    # Verify no has_key remains
    remaining = 0
    for py_file in sorted(root_path.rglob("*.py")):
        if "__pycache__" in str(py_file):
            continue
        content = py_file.read_text(encoding="utf-8", errors="replace")
        remaining += content.count(".has_key(")
    print(f"Remaining .has_key( occurrences: {remaining}")

    if remaining > 0:
        # Show the remaining ones for manual review
        for py_file in sorted(root_path.rglob("*.py")):
            if "__pycache__" in str(py_file):
                continue
            content = py_file.read_text(encoding="utf-8", errors="replace")
            for i, line in enumerate(content.splitlines(), 1):
                if ".has_key(" in line:
                    print(f"  UNFIXED {py_file}:{i}: {line.strip()}")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fix_has_key.py <directory>")
        sys.exit(1)
    main(sys.argv[1])
