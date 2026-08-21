#!/usr/bin/env python3

"""
Run an NDCH model in p4.

Input:
    alignment.nex
    clean_tree.nwk

The tree should contain four taxa:
    A, B, C, D

The AB clade receives a separate nucleotide composition.
"""

import json

# -----------------------------
# Input and output files
# -----------------------------

ALIGNMENT_FILE = "alignment.nex"
TREE_FILE = "clean_tree.nwk"
RESULT_FILE = "ndch_result.json"

# Starting value for the AB composition.
# The values must sum to 1.0.
AB_INITIAL_FREQS = [0.40, 0.10, 0.10, 0.40]


# -----------------------------
# Read alignment and tree
# -----------------------------

read(ALIGNMENT_FILE)
data = Data()

read(TREE_FILE)
t = var.trees[0]

# Attach the alignment to the tree.
t.data = data


# -----------------------------
# Find the AB internal node
# -----------------------------

def get_leaf_names(node):
    """Return all descendant tip names below a p4 node."""
    leaves = []

    def walk(n):
        if n.isLeaf:
            leaves.append(n.name)
            return

        child = n.leftChild
        while child is not None:
            walk(child)
            child = child.sibling

    walk(node)
    return sorted(leaves)


ab_parent = None

for node in t.iterNodesNoRoot():
    if get_leaf_names(node) == ["A", "B"]:
        ab_parent = node
        break

if ab_parent is None:
    raise RuntimeError("Could not find the internal node subtending A and B")

print("AB parent node number:", ab_parent.nodeNum)


# -----------------------------
# Define the NDCH model
# -----------------------------

# Composition 0:
# Equal nucleotide frequencies, fixed.
c_equal = t.newComp(
    partNum=0,
    free=0,
    spec="equal"
)

# Composition 1:
# Alternative nucleotide frequencies for the AB clade.
# free=1 means p4 will estimate these frequencies.
c_ab = t.newComp(
    partNum=0,
    free=1,
    spec="specified",
    val=AB_INITIAL_FREQS
)

# Put equal composition on the whole tree.
t.setModelComponentOnNode(
    c_equal,
    node=t.root,
    clade=1
)

# Override the composition on the AB clade.
t.setModelComponentOnNode(
    c_ab,
    node=ab_parent,
    clade=1
)


# -----------------------------
# Define the nucleotide
# substitution model
# -----------------------------

# HKY rate matrix.
# p4 uses spec="2p" for a HKY/K2P-style transition/transversion matrix.
# Here kappa is fixed at 2.0.
t.newRMatrix(
    partNum=0,
    free=0,
    spec="2p",
    val=2.0
)

# No gamma-rate variation.
t.setNGammaCat(
    partNum=0,
    nGammaCat=1
)

# No invariant-site category.
t.setPInvar(
    partNum=0,
    free=0,
    val=0.0
)


# -----------------------------
# Check and optimize
# -----------------------------

print("\nModel:")
t.model.dump()

print("\nOptimizing likelihood...")
t.optLogLike(
    verbose=1,
    method="BOBYQA",
    optBrLens=True
)

print("\nOptimized log likelihood:", t.logLike)


# -----------------------------
# Extract estimated parameters
# -----------------------------

result = {
    "logL": t.logLike,
    "compositions": {},
    "branch_lengths": {},
    "node_model_components": {}
}

# Save every composition in the model.
for comp in t.model.parts[0].comps:
    result["compositions"][str(comp.num)] = [
        float(x) for x in comp.val
    ]

# Save terminal and internal branch lengths.
for node in t.iterNodesNoRoot():
    if node.isLeaf:
        label = node.name
    elif node.nodeNum == ab_parent.nodeNum:
        label = "AB_internal"
    else:
        label = "internal_%d" % node.nodeNum

    result["branch_lengths"][label] = float(node.br.len)

# Record which composition is assigned to each node.
for node in t.iterNodes():
    if node.parts:
        result["node_model_components"][str(node.nodeNum)] = {
            "composition_number": int(node.parts[0].compNum),
            "is_leaf": bool(node.isLeaf),
            "name": node.name if node.isLeaf else None,
            "leaf_descendants": get_leaf_names(node)
        }


# -----------------------------
# Write results
# -----------------------------

with open(RESULT_FILE, "w") as f:
    json.dump(result, f, indent=2)

print("\nResults written to:", RESULT_FILE)

print("\nEstimated compositions:")
for comp_num, frequencies in result["compositions"].items():
    print("  comp", comp_num, frequencies)

print("\nEstimated branch lengths:")
for node_name, branch_length in result["branch_lengths"].items():
    print(" ", node_name, branch_length)
