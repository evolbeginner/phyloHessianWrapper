"""TNode/TTree structures and Newick-format tree parser."""

import numpy as np


class TNode:
    """Phylogenetic tree node with three branches (parent + two children)."""
    __slots__ = ('branch0', 'branch1', 'branch2', 'length0', 'length1', 'length2',
                 'param', 'comp_idx', 'rmat_idx', 'tip_no', 'sequence')

    def __init__(self):
        self.branch0 = None
        self.branch1 = None
        self.branch2 = None
        self.length0 = 0.0
        self.length1 = 0.0
        self.length2 = 0.0
        self.param = 1.0
        self.comp_idx = 0
        self.rmat_idx = 0
        self.tip_no = -1
        self.sequence = None


class TTree:
    """Phylogenetic tree with root, node list, and tip metadata."""
    __slots__ = ('root', 'node_list', 'num_tips', 'num_nodes', 'total_length',
                 'rooted', 'names', 'tips', 'partition_length', 'partition_rate')

    def __init__(self):
        self.root = None
        self.node_list = []
        self.num_tips = 0
        self.num_nodes = 0
        self.total_length = 0.0
        self.rooted = True
        self.names = {}
        self.tips = {}
        self.partition_length = -1
        self.partition_rate = 1.0


class ParseError(Exception):
    """Tree parsing error."""
    pass


def count_trees(text):
    """Count the number of trees (semicolons) in the input."""
    return text.count(';')


def parse_trees(text, translate_map=None):
    """Parse all trees from Newick-format text.

    Returns list of (partition_length, partition_rate, TTree) tuples.
    """
    parser = _NewickParser(text, translate_map or {})
    trees = []
    while parser._pos < len(text):
        parser._skip_space()
        while parser._peek() == ';':
            parser._advance()
            parser._skip_space()
        if parser._pos >= len(text):
            break
        trees.append(parser._parse_one_tree())
    return trees


class _NewickParser:
    """Recursive-descent Newick parser.

    Grammar:  tree = [partition_label] '(' subtree_list ')' [label] [':'] [metadata] ';'
              subtree_list = subtree {',' subtree}
              subtree = '(' subtree_list ')' [label] metadata | name metadata
              metadata = [':' float] ['[' annotation ']']
    """

    def __init__(self, text, translate_map):
        self._text = text
        self._pos = 0
        self._translate = translate_map

    # ---- character-level helpers ----

    def _peek(self):
        if self._pos >= len(self._text):
            return '\0'
        return self._text[self._pos]

    def _advance(self):
        if self._pos < len(self._text):
            self._pos += 1

    def _skip_space(self):
        while self._pos < len(self._text) and self._text[self._pos].isspace():
            self._pos += 1
        return self._peek()

    def _require(self, ch, msg):
        c = self._skip_space()
        if c != ch:
            raise ParseError(f"{msg}: expected '{ch}', got '{c}' at pos {self._pos}")
        self._advance()

    # ---- name reading ----

    def _parse_name(self):
        """Read an unquoted or single-quoted taxon name."""
        if self._peek() == "'":
            return self._parse_quoted_name()
        chars = []
        while self._pos < len(self._text):
            ch = self._text[self._pos]
            if ch in ':,) \t\n\r;[':
                break
            if ch not in ("'", '"'):
                chars.append(ch)
            self._pos += 1
        return ''.join(chars)

    def _parse_quoted_name(self):
        """Read a single-quoted name ('' escapes to ')."""
        self._advance()
        chars = []
        while self._pos < len(self._text):
            ch = self._text[self._pos]
            if ch == "'":
                if self._pos + 1 < len(self._text) and self._text[self._pos + 1] == "'":
                    chars.append("'")
                    self._pos += 2
                    continue
                self._advance()
                return ''.join(chars)
            chars.append(ch)
            self._pos += 1
        raise ParseError("unterminated quoted name")

    def _skip_label(self):
        """Consume an internal node label (bootstrap value, etc.) without storing it."""
        self._skip_space()
        while self._peek() not in (':', ',', ')', ';', '[', '\0'):
            self._pos += 1
            self._skip_space()

    # ---- numeric tokens ----

    def _read_number(self):
        """Read a floating-point token at the current position, or None."""
        self._skip_space()
        start = self._pos
        while self._pos < len(self._text) and self._text[self._pos] in '0123456789.-+eE':
            self._pos += 1
        if start == self._pos:
            return None
        token = self._text[start:self._pos]
        try:
            return float(token)
        except ValueError as exc:
            raise ParseError(f"invalid float '{token}' at pos {start}") from exc

    # ---- metadata (branch length + annotations) ----

    def _parse_metadata(self, node):
        """Read and apply :branch_length and [annotation] for a node.

        Returns the branch length (0.0 if absent).
        """
        blen = 0.0
        self._skip_space()
        if self._peek() == ':':
            self._advance()
            val = self._read_number()
            if val is None:
                raise ParseError("expected branch length after ':'")
            blen = val
            node.length0 = blen

        self._skip_space()
        if self._peek() == '[':
            self._advance()
            self._apply_annotation(node)

        return blen

    def _apply_annotation(self, node):
        """Parse [...] annotation content after a branch and apply to node."""
        buf = []
        while self._pos < len(self._text):
            ch = self._text[self._pos]
            if ch == ']':
                self._advance()
                break
            buf.append(ch)
            self._pos += 1
        content = ''.join(buf).strip()
        if not content:
            return
        for token in content.split(','):
            token = token.strip()
            if not token:
                continue
            if '=' in token:
                key, val = token.split('=', 1)
                key = key.strip().lower()
                try:
                    ival = int(val.strip())
                except ValueError:
                    raise ParseError(f"invalid annotation value: {token}")
                if key == 'f':
                    node.comp_idx = ival
                elif key == 'r':
                    node.rmat_idx = ival
                else:
                    raise ParseError(f"unknown annotation key: '{key}'")
            else:
                try:
                    node.param = float(token)
                except ValueError:
                    raise ParseError(f"invalid annotation: {token}")

    # ---- partition info ----

    def _parse_partition(self):
        """Read [length] or [length, rate] prefix before a tree."""
        self._require('[', "missing '[' for partition")
        while self._pos < len(self._text) and self._text[self._pos].isspace():
            self._pos += 1
        start = self._pos
        while self._pos < len(self._text) and self._text[self._pos].isdigit():
            self._pos += 1
        if start == self._pos:
            raise ParseError("expected partition length")
        p_len = int(self._text[start:self._pos])

        p_rate = 1.0
        self._skip_space()
        if self._peek() == ',':
            self._advance()
            val = self._read_number()
            if val is not None:
                p_rate = val
        self._skip_space()
        self._require(']', "missing ']' after partition information")
        return p_len, p_rate

    # ---- node construction ----

    def _make_terminal(self, tree):
        """Create a terminal (tip) node with its name."""
        node = TNode()
        name = self._parse_name()
        if self._translate and name.isdigit():
            mapped = self._translate.get(int(name))
            if mapped is not None:
                name = mapped
        if not name:
            raise ParseError(f"missing taxon name at pos {self._pos}")

        tip_no = tree.num_tips
        tree.names[tip_no] = name
        node.tip_no = tip_no
        tree.tips[tip_no] = node
        tree.num_tips += 1
        tree.node_list.append(node)
        tree.num_nodes += 1
        return node

    def _make_internal(self, tree):
        """Create an internal node and register it in the tree."""
        node = TNode()
        node.tip_no = -1
        tree.node_list.append(node)
        tree.num_nodes += 1
        return node

    # ---- recursive descent ----

    def _parse_subtree(self, tree, detect_polytomy):
        """Parse a single subtree: either ( internal_node ) or a terminal.

        Returns the parsed node.
        """
        ch = self._skip_space()
        if ch == '(':
            self._advance()
            node = self._parse_internal(tree, detect_polytomy)
            self._require(')', "missing ')'")
            self._skip_label()
        else:
            node = self._make_terminal(tree)

        # Attach branch metadata to the incoming branch (length0).
        self._parse_metadata(node)
        return node

    def _parse_internal(self, tree, detect_polytomy):
        """Parse the interior of a parenthesized group: child1, child2 [, ...].

        Reads exactly 2 children for rooted nodes; raises ParseError if
        detect_polytomy=True and a third child is encountered.
        """
        node = self._make_internal(tree)

        node.branch1 = self._parse_subtree(tree, detect_polytomy)
        node.branch1.branch0 = node
        node.length1 = node.branch1.length0

        self._require(',', "missing ',' between branches")

        node.branch2 = self._parse_subtree(tree, detect_polytomy)
        node.branch2.branch0 = node
        node.length2 = node.branch2.length0

        if detect_polytomy:
            self._skip_space()
            if self._peek() == ',':
                raise ParseError(
                    "tree contains polytomies; resolve with zero-length branches")

        return node

    def _parse_one_tree(self):
        """Parse one complete Newick tree.

        Returns (partition_length, partition_rate, TTree).
        """
        tree = TTree()
        p_len = -1
        p_rate = 1.0

        ch = self._skip_space()
        if ch == '[':
            p_len, p_rate = self._parse_partition()
            ch = self._skip_space()

        self._require('(', "expected '(' to start tree")

        # Root node subtree — read first two children without polytomy check
        tree.root = self._parse_internal(tree, detect_polytomy=False)

        # Check for third child → unrooted tree
        ch = self._skip_space()
        if ch == ',':
            self._advance()
            tree.rooted = False
            tree.root.branch0 = self._parse_subtree(tree, detect_polytomy=True)
            tree.root.branch0.branch0 = tree.root
            tree.root.length0 = tree.root.branch0.length0

        self._require(')', "missing ')' after tree")
        self._skip_label()
        self._parse_metadata(tree.root)  # root branch length + annotation

        self._skip_space()
        if self._peek() != ';':
            raise ParseError(f"expected ';' at end of tree at pos {self._pos}")
        self._advance()

        tree.total_length = _compute_tree_height(tree)
        tree.partition_length = p_len
        tree.partition_rate = p_rate

        return p_len, p_rate, tree


def _compute_tree_height(tree):
    """Compute the maximum root-to-tip path length (tree height).

    Uses an iterative stack-based depth-first search to avoid recursion limits.
    """
    max_height = 0.0
    stack = [(tree.root, 0.0)]
    while stack:
        node, acc = stack.pop()
        if node.tip_no >= 0:
            if acc > max_height:
                max_height = acc
            continue
        if node.branch1:
            stack.append((node.branch1, acc + node.branch1.length0))
        if node.branch2:
            stack.append((node.branch2, acc + node.branch2.length0))
        if not tree.rooted and node is tree.root and node.branch0:
            stack.append((node.branch0, acc + node.branch0.length0))
    return max_height


def allocate_sequences(tree, num_sites):
    """Allocate numpy arrays for all node sequences."""
    for node in tree.node_list:
        node.sequence = np.zeros(num_sites, dtype=np.int32)


def dispose_tree(tree):
    """Clear node list and sequences from a tree for reuse."""
    for node in tree.node_list:
        node.sequence = None
    tree.node_list.clear()
    tree.root = None
    tree.num_nodes = 0
    tree.num_tips = 0
    tree.total_length = 0.0
    tree.rooted = True
    tree.names = {}
    tree.tips = {}
    tree.partition_length = -1
    tree.partition_rate = 1.0
