"""NEXUS reader: extracts plain Newick trees, TRANSLATE map, and [&R]/[&U] rooted hints
from NEXUS-wrapped tree files (BEAST, MrBayes, etc.)."""

import re


def is_nexus(text):
    """Return True if text appears to be NEXUS format."""
    return bool(re.search(r'(?im)^\s*#NEXUS|^\s*BEGIN\s+\w+\s*;', text))


def _strip_bracket_content(text):
    """Remove all [...] enclosed content, handling nested brackets."""
    result = []
    depth = 0
    for ch in text:
        if ch == '[':
            depth += 1
        elif ch == ']':
            if depth > 0:
                depth -= 1
        elif depth == 0:
            result.append(ch)
    return ''.join(result)


def _parse_translate(block):
    """Parse TRANSLATE command: returns {int_id: name, ...}."""
    translate = {}
    m = re.search(r'(?is)TRANSLATE\s+(.+?);', block)
    if not m:
        return translate
    body = _strip_bracket_content(m.group(1))
    for item in body.split(','):
        item = item.strip()
        if not item:
            continue
        parts = item.split(maxsplit=1)
        if len(parts) == 2:
            try:
                translate[int(parts[0])] = parts[1]
            except ValueError:
                pass
    return translate


def _extract_tree_newicks(block):
    """Extract (clean_newick, rooted_flag) from a TREES block.
    rooted_flag: True ([&R]), False ([&U]), or None."""
    results = []
    clean = _strip_bracket_content(block)

    for m in re.finditer(r'(?i)(?:TREE|UTREE)\s+\*?\s*(\S+)\s*=', clean):
        tree_name = m.group(1)
        pos = m.end()

        paren_start = clean.find('(', pos)
        if paren_start == -1:
            continue

        depth = 0
        paren_end = None
        for i in range(paren_start, len(clean)):
            if clean[i] == '(':
                depth += 1
            elif clean[i] == ')':
                depth -= 1
                if depth == 0:
                    paren_end = i + 1
                    break

        if paren_end is None:
            continue

        semi_pos = clean.find(';', paren_end)
        if semi_pos == -1:
            continue

        newick = clean[paren_start:semi_pos + 1]

        rooted = _detect_rooted(block, tree_name)
        results.append((newick, rooted))

    return results


def _detect_rooted(original_block, tree_name):
    """Search between '=' and first '(' in the original block for [&R] or [&U]."""
    pattern = r'(?i)(?:TREE|UTREE)\s+\*?\s*' + re.escape(tree_name) + r'\s*='
    m = re.search(pattern, original_block)
    if not m:
        return None

    pos = m.end()
    depth = 0
    raw = []
    while pos < len(original_block):
        ch = original_block[pos]
        if ch == '[':
            depth += 1
        elif ch == ']':
            if depth > 0:
                depth -= 1
        elif ch == '(' and depth == 0:
            break
        raw.append(ch)
        pos += 1

    prefix = ''.join(raw)
    if '[&R]' in prefix:
        return True
    if '[&U]' in prefix:
        return False
    return None


def parse_nexus_trees(text):
    """Parse NEXUS text, returning (translate_map, [(newick_str, rooted_flag), ...]).
    rooted_flag: True, False, or None when no [&R]/[&U] marker present.
    Raises ValueError if no trees are found."""
    translate = {}
    all_trees = []

    blocks = re.findall(r'(?is)BEGIN\s+TREES\s*;(.*?)(?:END|ENDBLOCK)\s*;', text)

    for block in blocks:
        translate.update(_parse_translate(block))
        all_trees.extend(_extract_tree_newicks(block))

    if not all_trees:
        raise ValueError("No trees found in NEXUS file")
    return translate, all_trees
