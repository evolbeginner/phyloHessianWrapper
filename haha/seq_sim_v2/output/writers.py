"""Output formatters: PHYLIP, Relaxed PHYLIP, NEXUS, FASTA."""


def write_sequences(f, fmt, tree_set, partition_lengths, num_taxa, total_sites,
                    state_chars, tree_no=0, dataset_no=0, is_nuc=True):
    """Dispatch to the appropriate output writer."""
    if fmt == 'p':
        _write_phylip(f, tree_set, partition_lengths, num_taxa, total_sites,
                      state_chars)
    elif fmt == 'r':
        _write_relaxed_phylip(f, tree_set, partition_lengths, num_taxa,
                              total_sites, state_chars)
    elif fmt == 'n':
        _write_nexus(f, tree_set, partition_lengths, num_taxa, total_sites,
                     state_chars, tree_no, dataset_no, is_nuc)
    elif fmt == 'f':
        _write_fasta(f, tree_set, partition_lengths, num_taxa, total_sites,
                     state_chars)


def _write_phylip(f, tree_set, partition_lengths, num_taxa, total_sites,
                  state_chars):
    """Standard PHYLIP format: 10-char padded names, no leading space."""
    f.write(f"{num_taxa} {total_sites}\n")
    for tip_no in range(num_taxa):
        name = tree_set[0].names[tip_no][:10]
        f.write(f"{name:<10}")
        for k, tree in enumerate(tree_set):
            pt = partition_lengths[k]
            seq = tree.tips[tip_no].sequence[:pt]
            f.write(''.join(state_chars[int(s)] for s in seq))
        f.write('\n')


def _write_relaxed_phylip(f, tree_set, partition_lengths, num_taxa, total_sites,
                          state_chars):
    """Relaxed PHYLIP: names+space, no truncation, no leading space."""
    f.write(f"{num_taxa} {total_sites}\n")
    for tip_no in range(num_taxa):
        name = tree_set[0].names[tip_no]
        f.write(f"{name} ")
        for k, tree in enumerate(tree_set):
            pt = partition_lengths[k]
            seq = tree.tips[tip_no].sequence[:pt]
            f.write(''.join(state_chars[int(s)] for s in seq))
        f.write('\n')


def _write_nexus(f, tree_set, partition_lengths, num_taxa, total_sites,
                 state_chars, tree_no, dataset_no, is_nuc):
    """NEXUS format with DATA block (space indentation)."""
    if tree_no > 0 and dataset_no > 0:
        f.write(f"Begin DATA; [Tree {tree_no}, Dataset {dataset_no}]\n")
    elif tree_no > 0:
        f.write(f"Begin DATA; [Tree {tree_no}]\n")
    elif dataset_no > 0:
        f.write(f"Begin DATA; [Dataset {dataset_no}]\n")
    else:
        f.write("Begin DATA;\n")

    datatype = "DNA" if is_nuc else "PROTEIN"
    f.write(f"  Dimensions NTAX={num_taxa} NCHAR={total_sites};\n")
    f.write(f"  Format MISSING=? GAP=- DATATYPE={datatype};\n")
    f.write("  Matrix\n")

    max_len = max(len(tree_set[0].names[i]) for i in range(num_taxa))
    for tip_no in range(num_taxa):
        name = tree_set[0].names[tip_no]
        f.write(name)
        f.write(' ' * (max_len - len(name) + 1))
        for k, tree in enumerate(tree_set):
            pt = partition_lengths[k]
            seq = tree.tips[tip_no].sequence[:pt]
            f.write(''.join(state_chars[int(s)] for s in seq))
        f.write('\n')
    f.write("  ;\nEnd;\n\n")


def _write_fasta(f, tree_set, partition_lengths, num_taxa, total_sites,
                 state_chars):
    """FASTA format: >Name then 72-char wrapped sequence."""
    for tip_no in range(num_taxa):
        name = tree_set[0].names[tip_no]
        f.write(f">{name}\n")
        parts = []
        for k, tree in enumerate(tree_set):
            pt = partition_lengths[k]
            parts.append(''.join(state_chars[int(s)]
                          for s in tree.tips[tip_no].sequence[:pt]))
        seq = ''.join(parts)
        for i in range(0, len(seq), 72):
            f.write(seq[i:i + 72] + '\n')


def write_ancestral(f, fmt, tree, state_chars):
    """Write ancestral (internal node) sequences."""
    total_nodes = _count_nodes(tree.root)
    if not tree.rooted:
        total_nodes += _count_nodes(tree.root.branch0)
    total_sites = len(tree.root.sequence)

    if fmt == 'f':
        _write_anc_fasta(f, tree, state_chars)
    else:
        f.write(f"{total_nodes} {total_sites}\n")
        _write_anc_phylip(f, tree, state_chars, tree)


def _count_nodes(node):
    """Count all nodes in a subtree (follows branch1/branch2 only, never parent)."""
    if node is None:
        return 0
    return 1 + _count_nodes(node.branch1) + _count_nodes(node.branch2)


def _write_anc_fasta(f, tree, state_chars):
    """Write ancestral sequences in FASTA format (pre-order).
    Internal node numbering starts at num_tips + 1."""
    node_no = [tree.num_tips]
    _anc_fasta_node(f, tree.root, state_chars, tree, node_no)
    _anc_fasta_node(f, tree.root.branch1, state_chars, tree, node_no)
    _anc_fasta_node(f, tree.root.branch2, state_chars, tree, node_no)
    if not tree.rooted:
        _anc_fasta_node(f, tree.root.branch0, state_chars, tree, node_no)


def _anc_fasta_node(f, node, state_chars, tree, node_no):
    """Recursive FASTA ancestral writer."""
    if node is None:
        return
    if node.tip_no == -1:
        node_no[0] += 1
        label = str(node_no[0])
    else:
        label = tree.names[node.tip_no]
    seq = ''.join(state_chars[int(s)] for s in node.sequence)
    f.write(f">{label}\n")
    for i in range(0, len(seq), 72):
        f.write(seq[i:i + 72] + '\n')
    if node.tip_no == -1:
        _anc_fasta_node(f, node.branch1, state_chars, tree, node_no)
        _anc_fasta_node(f, node.branch2, state_chars, tree, node_no)


def _write_anc_phylip(f, tree, state_chars, root_tree):
    """Write ancestral sequences in PHYLIP format.
    Internal node numbering starts at num_tips + 1."""
    node_no = [root_tree.num_tips]
    _anc_phy_node(f, tree.root, state_chars, root_tree, node_no)
    _anc_phy_node(f, tree.root.branch1, state_chars, root_tree, node_no)
    _anc_phy_node(f, tree.root.branch2, state_chars, root_tree, node_no)
    if not tree.rooted:
        _anc_phy_node(f, tree.root.branch0, state_chars, root_tree, node_no)


def _anc_phy_node(f, node, state_chars, tree, node_no):
    """Recursive PHYLIP ancestral writer."""
    if node is None:
        return
    if node.tip_no == -1:
        node_no[0] += 1
        label = str(node_no[0])
    else:
        label = tree.names[node.tip_no]
    seq = ''.join(state_chars[int(s)] for s in node.sequence)
    f.write(f"{label:<10}{seq}\n")
    if node.tip_no == -1:
        _anc_phy_node(f, node.branch1, state_chars, tree, node_no)
        _anc_phy_node(f, node.branch2, state_chars, tree, node_no)


def write_rates(f, site_rates, num_sites):
    """Write site-specific relative substitution rates."""
    f.write("Relative rates for each site:\n")
    for i in range(num_sites):
        f.write(f"{i + 1}\t{site_rates[i]:.6f}\n")
    f.write('\n')
