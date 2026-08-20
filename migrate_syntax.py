#!/usr/bin/env python3
"""
Python 2.7 -> 3.13 syntax migration — comprehensive script.

Handles:
1. print 'x' -> print('x')  (skips multi-line expressions with unbalanced parens/quotes)
2. print >> file, 'x' -> print('x', file=file)
3. print '''...''' -> print('''...''')  (single and multi-line triple-quoted)
4. except E, e: -> except E as e:
5. apply(func, args) -> func(*args)
6. raise ExcType, msg -> raise ExcType(msg)
7. raise exc_info[0], exc_info[1], exc_info[2] -> raise exc_info[1]
"""
import re
import sys
from pathlib import Path


def has_unbalanced(s: str) -> bool:
    """Check if line has unbalanced parens or quotes (expression continues)."""
    depth = 0
    in_sq = in_dq = False
    i = 0
    while i < len(s):
        c = s[i]
        if c == '\\' and i + 1 < len(s):
            i += 2
            continue
        if in_sq:
            if c == "'":
                in_sq = False
        elif in_dq:
            if c == '"':
                in_dq = False
        else:
            if c == "'":
                in_sq = True
            elif c == '"':
                in_dq = True
            elif c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
        i += 1
    return depth != 0 or in_sq or in_dq


def fix_line(line: str) -> str:
    """Apply fixes to a single line. Returns the fixed line (with original newline)."""
    newline = ''
    if line.endswith('\n'):
        newline = '\n'
        cp = line[:-1]
    else:
        cp = line

    stripped = cp.lstrip()
    indent = cp[:len(cp) - len(cp.lstrip())]

    # --- Skip comments, empty, already-correct print() ---
    if stripped.startswith('#') or not stripped:
        return line
    if stripped.startswith("print(") or stripped.startswith("print ("):
        return line

    # --- Skip lines with unbalanced parens (multi-line expressions) ---
    # EXCEPT print >> (which should be self-contained)
    if not stripped.startswith("print >>") and has_unbalanced(cp):
        return line

    # --- print >> file, ... ---
    m = re.match(r'^print\s*>>\s*(\w+)\s*,\s*(\S.*)$', stripped)
    if m:
        target, rest = m.groups()
        return f'{indent}print({rest}, file={target}){newline}'

    # --- bare print ---
    if stripped == 'print':
        return f'{indent}print(){newline}'

    # --- print 'x' / print a, b, c ---
    m = re.match(r'^print\s+(.+)$', stripped)
    if m:
        rest = m.group(1).rstrip()
        if rest:
            return f'{indent}print({rest}){newline}'
        return f'{indent}print(){newline}'

    # --- except E, e: ---
    m = re.match(r'^except\s+([\w.]+)\s*,\s*(\w+)\s*:', stripped)
    if m:
        exc, var = m.groups()
        rest = stripped[m.end():]
        return f'{indent}except {exc} as {var}:{rest}{newline}'

    # --- raise exc_info[0], exc_info[1], exc_info[2] ---
    if re.match(r'^raise\s+exc_info\[0\],\s*exc_info\[1\],\s*exc_info\[2\]$', stripped):
        return f'{indent}raise exc_info[1]{newline}'

    # --- raise ExcType, msg (2-arg) ---
    m = re.match(r'^raise\s+(\w+)\s*,\s*(.+)$', stripped)
    if m and not stripped.startswith("raise ") or (m and '(' not in m.group(1)):
        exc = m.group(1)
        rest = m.group(2).rstrip()
        if not rest.startswith('from '):
            return f'{indent}raise {exc}({rest}){newline}'

    # --- apply(func, args) ---
    m = re.match(r'^apply\s*\(\s*(\w+)\s*,\s*(.+)\)$', stripped)
    if m:
        func, args = m.groups()
        return f'{indent}{func}(*{args.rstrip()}){newline}'

    return line


def fix_triple_quoted_print(content: str) -> str:
    """Fix print triple-quoted strings (single and multi-line)."""
    # Single-line triple-quoted
    content = re.sub(
        r'^(\s*)print\s+("""[^"].*?""")\s*$',
        r'\1print(\2)',
        content, flags=re.MULTILINE)
    content = re.sub(
        r'^(\s*)print\s+(\x27\x27\x27[^\x27]*?\x27\x27\x27)\s*$',
        r'\1print(\2)',
        content, flags=re.MULTILINE)

    # Multi-line triple-quoted: print ''' followed by content until '''
    def multi_repl(m):
        indent = m.group(1)
        text = m.group(2)
        return f'{indent}print({text})'

    content = re.sub(
        r'^(\s*)print\s+(\"\"\".*?\"\"\")\s*$',
        multi_repl,
        content, flags=re.DOTALL | re.MULTILINE)
    content = re.sub(
        r'^(\s*)print\s+(\x27\x27\x27.*?\x27\x27\x27)\s*$',
        multi_repl,
        content, flags=re.DOTALL | re.MULTILINE)

    return content


def migrate_file(filepath: Path) -> bool:
    try:
        content = filepath.read_text(encoding='utf-8')
    except (UnicodeDecodeError, OSError):
        return False

    original = content

    # Step 1: fix triple-quoted print (multi-line safe)
    content = fix_triple_quoted_print(content)

    # Step 2: fix single-line print/except/raise/apply (skip unbalanced)
    lines = content.splitlines(keepends=True)
    lines = [fix_line(l) for l in lines]
    content = ''.join(lines)

    if content != original:
        filepath.write_text(content, encoding='utf-8')
        return True
    return False


def main():
    if len(sys.argv) < 2:
        print("Usage: migrate_syntax.py <root_dir>")
        sys.exit(1)

    root = Path(sys.argv[1])
    modified = 0
    total = 0

    for filepath in sorted(root.rglob("*.py")):
        if any(part in ('.venv', 'build', '__pycache__') for part in filepath.parts):
            continue
        total += 1
        if migrate_file(filepath):
            modified += 1

    print(f"Processed {total} files, modified {modified}")


if __name__ == "__main__":
    main()
