#!/usr/bin/env python3

"""
Unrooted NDCH analysis with command-line branch assignment.

Examples
--------

Optimize two NDCH composition vectors with equal exchangeabilities:

    python run_unrooted_ndch.py \
        -m JC+G \
        --comp_clade "t1,t12" \
        --pre JC \
        --alignment example-dna/sim/alignment/combined.fas \
        --tree example-dna/ref.tre

Select all branches on a path:

    python run_unrooted_ndch.py \
        -m GTR+G \
        --comp 't1,t2|t5'

Select multiple paths:

    python run_unrooted_ndch.py \
        -m GTR+G \
        --comp 't1,t2|t5;t3,t4|t5'

Select all branches below an LCA:

    python run_unrooted_ndch.py \
        -m GTR+G \
        --comp_clade 't1,t9'

Fix all composition vectors to equal frequencies:

    python run_unrooted_ndch.py \
        -m JC+G \
        --comp_clade 't1,t12' \
        --equal-comp-freq

Use p4 directly:

    p4 run_unrooted_ndch.py -- \
        -m JC+G \
        --comp_clade "t1,t12" \
        --pre JC \
        --alignment example-dna/sim/alignment/combined.fas \
        --tree example-dna/ref.tre


NDCH interpretation
-------------------

Component 1 is assigned to all background branches.

Component 2 is assigned to branches selected by --comp and/or
--comp_clade.

Unless --equal-comp-freq is supplied, both composition vectors are
free and are optimized by p4.

Important model-name interpretation
-----------------------------------

When a model such as JC, K2P, or another equal-frequency IQ-TREE model
is used without --equal-comp-freq, only its exchangeability/rate
constraints are retained. The NDCH composition vectors remain free.

For example:

    -m JC+G

means:

    * equal exchangeabilities;
    * gamma-distributed rates;
    * independently optimized NDCH composition vectors.

This is intentionally not strict homogeneous JC, because strict JC
requires all nucleotide frequencies to remain equal.

Use:

    --equal-comp-freq

if strict equal nucleotide frequencies are desired.


--comp semantics
----------------

Each side of "|" is a comma-separated taxon group. The LCA of each
group is found, and every branch on the unique path connecting the two
LCA nodes is assigned component 2.

For example, for:

    (((t1,t2),(t3,t4)),t5);

the specification:

    --comp 't1,t2|t5'

selects every branch on the path from LCA(t1,t2) to tip t5.

Multiple specifications are separated by semicolons:

    --comp 't1,t2|t5;t3,t4|t5'


--comp_clade semantics
----------------------

Each specification is a comma-separated taxon group. The LCA of that
group is found under p4's current root orientation, and every
descendant branch inside that LCA clade is assigned component 2.

Multiple specifications are separated by semicolons:

    --comp_clade 't1,t12;t4,t7'
"""

from p4 import *

import argparse
import json
import math
import sys

from iqtree_models import configure_model, model_help_names, parse_model


# ============================================================
# Command-line parsing
# ============================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Run an unrooted NDCH analysis with one background "
            "composition and one composition assigned to selected "
            "branches."
        )
    )

    parser.add_argument(
        "-m",
        "--model",
        default="HKY+G",
        help=(
            "Substitution model in IQ-TREE notation. "
            "Default: HKY+G. Available base models: %s. "
            "Modifiers include +I, +G, and +G<number>. "
            "Unless --equal-comp-freq is used, NDCH composition "
            "vectors remain free even for equal-frequency base models."
            % model_help_names()
        )
    )

    parser.add_argument(
        "-s",
        "--alignment",
        default="alignment.nex",
        help=(
            "Input alignment readable by p4. "
            "Default: alignment.nex."
        )
    )

    parser.add_argument(
        "-t",
        "--tree",
        default="tree.nwk",
        help=(
            "Input Newick or Nexus tree readable by p4. "
            "Default: tree.nwk."
        )
    )

    parser.add_argument(
        "--pre",
        default="iqtree",
        help=(
            "Output prefix. Default: iqtree."
        )
    )

    parser.add_argument(
        "--comp",
        default="",
        help=(
            "Semicolon-separated path specifications assigned to "
            "component 2. Each side of '|' is a comma-separated "
            "taxon group. Example: 't1,t2|t5'."
        )
    )

    parser.add_argument(
        "--comp_clade",
        "--comp-clade",
        default="",
        help=(
            "Semicolon-separated unrooted clades assigned to "
            "component 2. Example: 't1,t12;t4,t7'."
        )
    )

    parser.add_argument(
        "--equal-comp-freq",
        action="store_true",
        help=(
            "Fix all composition vectors to "
            "A=C=G=T=0.25. By default they are optimized."
        )
    )

    parser.add_argument(
        "--optimizer",
        choices=[
            "BOBYQA",
            "newtAndBrentPowell",
            "newtAndBOBYQA",
            "allBrentPowell"
        ],
        default="BOBYQA",
        help=(
            "p4 likelihood optimizer. Default: BOBYQA."
        )
    )

    parser.add_argument(
        "--refine",
        action="store_true",
        help=(
            "After the requested optimizer, run a second optimization "
            "with newtAndBrentPowell. By default, only the requested "
            "optimizer is run."
        )
    )

    p4_args = getattr(var, "argvAfterDoubleDash", None)

    if p4_args is not None:
        command_line_args = list(p4_args)
    else:
        command_line_args = sys.argv[1:]

    return parser.parse_args(command_line_args)


args = parse_arguments()


ALIGNMENT_FILE = args.alignment
TREE_FILE = args.tree

TREE_OUTPUT_FILE = args.pre + ".treefile"
COMPONENT_TREE_OUTPUT_FILE = args.pre + ".comp.treefile"
RESULT_FILE = args.pre + ".result.json"
REPORT_FILE = args.pre + ".ndch"


# ============================================================
# Read alignment and tree
# ============================================================

read(ALIGNMENT_FILE)
d = Data()

read(TREE_FILE)
t = var.trees[0]

t.data = d

expected_taxa = list(d.taxNames)

if sorted(t.taxNames) != sorted(expected_taxa):
    raise RuntimeError(
        "Tree taxa do not match alignment taxa. "
        "Alignment taxa: %s; tree taxa: %s"
        % (
            sorted(expected_taxa),
            sorted(t.taxNames)
        )
    )


# ============================================================
# General helper functions
# ============================================================

def descendant_taxa(node):
    """
    Return the sorted names of all tips descended from node.
    """
    names = []

    def walk(current):
        if current.isLeaf:
            names.append(current.name)
            return

        child = current.leftChild

        while child is not None:
            walk(child)
            child = child.sibling

    walk(node)

    return sorted(names)


def node_label(node):
    """
    Return a readable node label.
    """
    if node is None:
        return "none"

    if node.isLeaf:
        return node.name

    return "node_%d" % int(node.nodeNum)


def branch_label(node):
    """
    Return a label for the branch represented by its child node.
    """
    return "%s -> %s" % (
        node_label(node.parent),
        node_label(node)
    )


def display_comp_num(internal_comp_num):
    """
    Convert p4's zero-based composition number to a one-based label.
    """
    return int(internal_comp_num) + 1


def iter_children(node):
    """
    Yield the immediate children of node.
    """
    child = node.leftChild

    while child is not None:
        yield child
        child = child.sibling


def descendant_branch_nodes(node):
    """
    Return all branch-representing nodes below node.
    """
    result = []

    def walk(current):
        for child in iter_children(current):
            result.append(child)
            walk(child)

    walk(node)

    return result


def find_tip(tree, taxon_name):
    """
    Find a tree tip by name.
    """
    for node in tree.iterNodesNoRoot():
        if node.isLeaf and node.name == taxon_name:
            return node

    raise RuntimeError(
        "Could not find tip '%s'." % taxon_name
    )


def ancestor_path(node):
    """
    Return node followed by all its ancestors through the root.
    """
    result = []
    current = node

    while current is not None:
        result.append(current)
        current = current.parent

    return result


def find_lca_node(tree, taxon_names):
    """
    Return the LCA of taxon_names under p4's current root orientation.
    """
    if not taxon_names:
        raise ValueError(
            "At least one taxon is required to calculate an LCA."
        )

    tips = [
        find_tip(tree, taxon_name)
        for taxon_name in taxon_names
    ]

    paths = [
        ancestor_path(tip)
        for tip in tips
    ]

    common_node_ids = set(
        id(node)
        for node in paths[0]
    )

    for path in paths[1:]:
        common_node_ids.intersection_update(
            id(node)
            for node in path
        )

    for node in paths[0]:
        if id(node) in common_node_ids:
            return node

    raise RuntimeError(
        "Could not calculate LCA(%s)."
        % ",".join(sorted(taxon_names))
    )


def branches_between_nodes(first_node, second_node):
    """
    Return all branch-representing child nodes on the unique path
    between first_node and second_node.
    """
    first_path = ancestor_path(first_node)
    second_path = ancestor_path(second_node)

    second_ids = set(
        id(node)
        for node in second_path
    )

    common_ancestor = None

    for node in first_path:
        if id(node) in second_ids:
            common_ancestor = node
            break

    if common_ancestor is None:
        raise RuntimeError(
            "No common ancestor found for %s and %s."
            % (
                node_label(first_node),
                node_label(second_node)
            )
        )

    branches = []

    current = first_node

    while current is not common_ancestor:
        branches.append(current)
        current = current.parent

    downward_branches = []
    current = second_node

    while current is not common_ancestor:
        downward_branches.append(current)
        current = current.parent

    branches.extend(reversed(downward_branches))

    return branches


def split_branch_spec(specification):
    """
    Split a path specification into left and right groups.
    """
    specification = specification.strip()

    if not specification:
        raise ValueError(
            "An empty path specification was supplied."
        )

    pieces = specification.split("|")

    if len(pieces) != 2:
        raise ValueError(
            "Path specification must contain exactly one '|': %s"
            % specification
        )

    left = pieces[0].strip()
    right = pieces[1].strip()

    if not left or not right:
        raise ValueError(
            "Both sides of the path specification are required: %s"
            % specification
        )

    return left, right


def expand_taxon_group(group, valid_taxa):
    """
    Convert a taxon-group string into a sorted list of taxa.
    """
    group = group.strip()

    if not group:
        raise ValueError(
            "An empty taxon group was supplied."
        )

    if "," in group:
        raw_taxa = group.split(",")

        if any(not item.strip() for item in raw_taxa):
            raise ValueError(
                "Taxon group '%s' contains an empty taxon name."
                % group
            )

        taxa = [
            item.strip()
            for item in raw_taxa
        ]

    elif group in valid_taxa:
        taxa = [group]

    else:
        # Backward-compatible syntax for single-character taxon names,
        # such as AB meaning A,B.
        taxa = list(group)

    unknown = [
        taxon
        for taxon in taxa
        if taxon not in valid_taxa
    ]

    if unknown:
        raise ValueError(
            "Unknown taxon name(s) in '%s': %s"
            % (
                group,
                ", ".join(unknown)
            )
        )

    if len(set(taxa)) != len(taxa):
        raise ValueError(
            "Taxon repeated in group '%s'." % group
        )

    return sorted(taxa)


def split_requested_specs(text):
    """
    Split semicolon-separated path specifications.
    """
    text = text.strip()

    if not text:
        return []

    if ";" in text:
        return [
            item.strip()
            for item in text.split(";")
            if item.strip()
        ]

    if text.count("|") == 1:
        return [text]

    # Legacy separator for specifications with single-character taxa.
    return [
        item.strip()
        for item in text.split(",")
        if item.strip()
    ]


def split_requested_clades(text):
    """
    Split semicolon-separated clade specifications.
    """
    text = text.strip()

    if not text:
        return []

    return [
        item.strip()
        for item in text.split(";")
        if item.strip()
    ]


# ============================================================
# Branch-selection functions
# ============================================================

def find_selected_branches(tree, specification):
    """
    Resolve one --comp path specification.
    """
    left_text, right_text = split_branch_spec(specification)

    valid_taxa = set(tree.taxNames)

    left_taxa = expand_taxon_group(
        left_text,
        valid_taxa
    )

    right_taxa = expand_taxon_group(
        right_text,
        valid_taxa
    )

    left_lca = find_lca_node(
        tree,
        left_taxa
    )

    right_lca = find_lca_node(
        tree,
        right_taxa
    )

    if left_lca is right_lca:
        raise ValueError(
            "Path '%s' resolves both groups to the same LCA (%s). "
            "The resulting path contains no branches."
            % (
                specification,
                node_label(left_lca)
            )
        )

    connecting_branches = branches_between_nodes(
        left_lca,
        right_lca
    )

    if not connecting_branches:
        raise RuntimeError(
            "Path '%s' did not produce any branches."
            % specification
        )

    print(
        "\nResolved path specification: %s"
        % specification
    )

    print(
        "  Left:  LCA(%s) = %s"
        % (
            left_text,
            node_label(left_lca)
        )
    )

    print(
        "  Right: LCA(%s) = %s"
        % (
            right_text,
            node_label(right_lca)
        )
    )

    print(
        "  Selected %d branch%s:"
        % (
            len(connecting_branches),
            "" if len(connecting_branches) == 1 else "es"
        )
    )

    for index, node in enumerate(
        connecting_branches,
        start=1
    ):
        print(
            "    %d. %s; descendants=%s"
            % (
                index,
                branch_label(node),
                ",".join(descendant_taxa(node))
            )
        )

    return {
        "specification": specification,
        "left_text": left_text,
        "right_text": right_text,
        "left_taxa": left_taxa,
        "right_taxa": right_taxa,
        "left_lca": left_lca,
        "right_lca": right_lca,
        "branches": connecting_branches
    }


def find_clade_branches(tree, clade_text):
    """
    Resolve one --comp_clade specification.

    The requested taxa define an LCA under p4's current root
    orientation. Every descendant branch inside that LCA clade is
    assigned component 2.
    """
    valid_taxa = set(tree.taxNames)

    clade_taxa = expand_taxon_group(
        clade_text,
        valid_taxa
    )

    if len(clade_taxa) < 2:
        raise ValueError(
            "LCA clade '%s' must contain at least two taxa."
            % clade_text
        )

    lca_node = find_lca_node(
        tree,
        clade_taxa
    )

    selected = descendant_branch_nodes(lca_node)

    if not selected:
        raise ValueError(
            "LCA clade '%s' contains no descendant branches."
            % clade_text
        )

    print(
        "\nResolved clade specification: %s"
        % clade_text
    )

    print(
        "  Rooted LCA: %s"
        % node_label(lca_node)
    )

    print(
        "  Rooted LCA descendants: %s"
        % ",".join(descendant_taxa(lca_node))
    )

    print(
        "  Selected %d branch%s:"
        % (
            len(selected),
            "" if len(selected) == 1 else "es"
        )
    )

    for index, node in enumerate(
        selected,
        start=1
    ):
        print(
            "    %d. %s; descendants=%s"
            % (
                index,
                branch_label(node),
                ",".join(descendant_taxa(node))
            )
        )

    return {
        "specification": clade_text,
        "taxa": clade_taxa,
        "lca": lca_node,
        "branches": selected
    }


def register_selected_branch(
    records,
    node,
    source_type,
    source_specification
):
    """
    Register a selected branch without duplicating overlapping paths.
    """
    node_number = int(node.nodeNum)

    if node_number not in records:
        records[node_number] = {
            "node": node,
            "sources": []
        }

    source = {
        "type": source_type,
        "specification": source_specification
    }

    if source not in records[node_number]["sources"]:
        records[node_number]["sources"].append(source)


def format_selection_sources(sources):
    """
    Return a readable branch-selection source string.
    """
    return "; ".join(
        "%s '%s'" % (
            source["type"],
            source["specification"]
        )
        for source in sources
    )


# ============================================================
# Parse and resolve requested branches
# ============================================================

requested_specs = split_requested_specs(args.comp)
requested_clades = split_requested_clades(args.comp_clade)


if len(set(requested_specs)) != len(requested_specs):
    raise ValueError(
        "The same path specification was supplied more than once."
    )

if len(set(requested_clades)) != len(requested_clades):
    raise ValueError(
        "The same clade specification was supplied more than once."
    )


print("Requested component 2 paths:")

if requested_specs:
    for specification in requested_specs:
        print("  %s" % specification)
else:
    print("  none")


print("Requested component 2 clades:")

if requested_clades:
    for specification in requested_clades:
        print("  %s" % specification)
else:
    print("  none")


resolved_paths = []
resolved_clades = []
selected_branch_records = {}


for specification in requested_specs:
    resolved_path = find_selected_branches(
        t,
        specification
    )

    resolved_paths.append(resolved_path)

    for selected_node in resolved_path["branches"]:
        register_selected_branch(
            selected_branch_records,
            selected_node,
            "path",
            specification
        )


for specification in requested_clades:
    resolved_clade = find_clade_branches(
        t,
        specification
    )

    resolved_clades.append(resolved_clade)

    for selected_node in resolved_clade["branches"]:
        register_selected_branch(
            selected_branch_records,
            selected_node,
            "clade",
            specification
        )


selected_branch_records = dict(
    sorted(
        selected_branch_records.items(),
        key=lambda item: item[0]
    )
)


print(
    "\nTotal unique branches assigned to component 2: %d"
    % len(selected_branch_records)
)


# ============================================================
# Parse the substitution model
# ============================================================

try:
    model_definition = parse_model(args.model)

except ValueError as error:
    raise SystemExit(
        "error: %s" % error
    )


model_name = model_definition["name"]
base_model = model_definition["base_model"]


# ============================================================
# Create NDCH composition components
# ============================================================

# The values supplied here are starting values when free=1.
#
# configure_model() is called later and may apply the equal-frequency
# constraints associated with models such as JC. Therefore the NDCH
# free/fixed status and starting values are explicitly restored after
# configure_model().

composition_free = 0 if args.equal_comp_freq else 1

background_start = [
    0.25,
    0.25,
    0.25,
    0.25
]

selected_start = [
    0.40,
    0.10,
    0.10,
    0.40
]


comp_0 = t.newComp(
    partNum=0,
    free=composition_free,
    spec="specified",
    val=background_start
)


comp_1 = None

if selected_branch_records:
    comp_1 = t.newComp(
        partNum=0,
        free=composition_free,
        spec="specified",
        val=(
            background_start
            if args.equal_comp_freq
            else selected_start
        )
    )


# ============================================================
# Assign composition components to branches
# ============================================================

# Set the background composition over the complete tree.
t.setModelComponentOnNode(
    comp_0,
    node=t.root,
    clade=1
)


# Override selected branches with component 2.
for branch_record in selected_branch_records.values():
    selected_node = branch_record["node"]

    print(
        "Assigning component 2 to %s; descendants=%s; selected_by=%s"
        % (
            branch_label(selected_node),
            ",".join(descendant_taxa(selected_node)),
            format_selection_sources(branch_record["sources"])
        )
    )

    t.setModelComponentOnNode(
        comp_1,
        node=selected_node,
        clade=0
    )


# ============================================================
# Configure rate matrix and rate heterogeneity
# ============================================================

configure_model(
    t,
    model_definition
)


# ============================================================
# Restore and verify NDCH composition settings
# ============================================================

def set_comp_values(comp, values):
    """
    Set composition values without replacing p4's internal val object.
    """
    if len(comp.val) != len(values):
        raise RuntimeError(
            "Composition dimension mismatch: p4 has %d values, "
            "but %d values were supplied."
            % (
                len(comp.val),
                len(values)
            )
        )

    for index, value in enumerate(values):
        comp.val[index] = float(value)


def restore_ndch_composition_settings():
    """
    Restore NDCH composition settings after configure_model().

    This is necessary because an IQ-TREE model such as JC may imply
    equal fixed frequencies. In this script, model names control the
    exchangeability constraints, while NDCH composition vectors remain
    independently optimizable unless --equal-comp-freq is supplied.
    """
    expected_free = 0 if args.equal_comp_freq else 1

    comp_0.free = expected_free
    set_comp_values(
        comp_0,
        background_start
    )

    if comp_1 is not None:
        comp_1.free = expected_free

        set_comp_values(
            comp_1,
            (
                background_start
                if args.equal_comp_freq
                else selected_start
            )
        )


restore_ndch_composition_settings()


starting_compositions = {
    int(comp.num): [
        float(value)
        for value in comp.val
    ]
    for comp in t.model.parts[0].comps
}


print("\nComposition settings before modelSanityCheck():")

for comp in t.model.parts[0].comps:
    print(
        "  comp %d: free=%d, values=%s"
        % (
            display_comp_num(comp.num),
            int(comp.free),
            [
                float(value)
                for value in comp.val
            ]
        )
    )


# ============================================================
# Check model and verify free status
# ============================================================

# resetEmpiricalComps=False prevents p4 from resetting composition
# vectors during model checking. These comps use spec='specified',
# but this also makes the intended behavior explicit.
t.modelSanityCheck(
    resetEmpiricalComps=False
)


expected_comp_free = 0 if args.equal_comp_freq else 1

for comp in t.model.parts[0].comps:
    if int(comp.free) != expected_comp_free:
        raise RuntimeError(
            "Composition component %d has free=%d after "
            "modelSanityCheck(), but free=%d was requested."
            % (
                display_comp_num(comp.num),
                int(comp.free),
                expected_comp_free
            )
        )


print("\nComposition settings after modelSanityCheck():")

for comp in t.model.parts[0].comps:
    print(
        "  comp %d: free=%d, values=%s"
        % (
            display_comp_num(comp.num),
            int(comp.free),
            [
                float(value)
                for value in comp.val
            ]
        )
    )


print(
    "\nNumber of free model parameters reported by p4: %d"
    % int(t.model.nFreePrams)
)


if args.equal_comp_freq:
    print(
        "Composition frequencies are fixed because "
        "--equal-comp-freq was supplied."
    )
else:
    print(
        "Composition frequencies are free and will be optimized."
    )


print("\nTree with model assignments:")
t.draw(model=1)


# ============================================================
# Optimize likelihood
# ============================================================

print(
    "\nOptimizing likelihood with %s ..."
    % args.optimizer
)

t.optLogLike(
    verbose=1,
    method=args.optimizer
)

refine_likelihood = args.refine

if refine_likelihood:
    print(
        "\nRefining likelihood with newtAndBrentPowell ..."
    )

    t.optLogLike(
        verbose=1,
        method="newtAndBrentPowell"
    )


print(
    "\nOptimized log likelihood: %.10f"
    % float(t.logLike)
)


# ============================================================
# Examine optimized composition vectors
# ============================================================

optimized_compositions = {
    int(comp.num): [
        float(value)
        for value in comp.val
    ]
    for comp in t.model.parts[0].comps
}


def maximum_absolute_change(start_values, end_values):
    """
    Return the largest absolute element-wise change.
    """
    return max(
        abs(end - start)
        for start, end in zip(start_values, end_values)
    )


print("\nStarting and optimized composition vectors:")

for comp in t.model.parts[0].comps:
    comp_number = int(comp.num)

    start_values = starting_compositions[comp_number]
    end_values = optimized_compositions[comp_number]

    max_change = maximum_absolute_change(
        start_values,
        end_values
    )

    print(
        "\n  Component %d"
        % display_comp_num(comp_number)
    )

    print(
        "    free:      %d"
        % int(comp.free)
    )

    print(
        "    starting:  %s"
        % [
            float(value)
            for value in start_values
        ]
    )

    print(
        "    optimized: %s"
        % [
            float(value)
            for value in end_values
        ]
    )

    print(
        "    maximum absolute change: %.12g"
        % max_change
    )

    if not args.equal_comp_freq and max_change < 1.0e-10:
        print(
            "    WARNING: this composition did not move detectably "
            "from its starting values."
        )


# ============================================================
# Substitution model and Q matrices
# ============================================================

part = t.model.parts[0]
r_matrix = part.rMatrices[0]
gdasrv = part.gdasrvs[0] if part.gdasrvs else None


print("\nSubstitution model:")
print("  requested model: %s" % model_name)
print("  base model:      %s" % base_model)

print(
    "  rate matrix:     %d"
    % int(r_matrix.num)
)

print(
    "  rate spec:       %s"
    % r_matrix.spec
)

print(
    "  rate free:       %d"
    % int(r_matrix.free)
)

print(
    "  rate values:     %s"
    % [
        float(value)
        for value in r_matrix.val
    ]
)

print(
    "  IQ-TREE rate constraint: %s"
    % model_definition["rate_code"]
)

print(
    "  p4 rate handling:        %s"
    % model_definition["rate_mode"]
)

print(
    "  invariant proportion: %.10g; free=%d"
    % (
        float(part.pInvar.val),
        int(part.pInvar.free)
    )
)


if not args.equal_comp_freq:
    print(
        "  NDCH frequencies: independently optimized"
    )

    if model_definition["equal_frequencies_in_iqtree"]:
        print(
            "  Note: the requested base model normally has equal "
            "frequencies, but this NDCH analysis overrides that "
            "frequency constraint."
        )

else:
    print(
        "  NDCH frequencies: fixed equal frequencies"
    )


if gdasrv is not None:
    print(
        "  gamma categories: %d"
        % int(part.nGammaCat)
    )

    print(
        "  gamma alpha:      %.10g; free=%d"
        % (
            float(gdasrv.val[0]),
            int(gdasrv.free)
        )
    )

    print(
        "  gamma frequencies: %s"
        % [
            float(value)
            for value in gdasrv.freqs
        ]
    )

    print(
        "  gamma rates:       %s"
        % [
            float(value)
            for value in gdasrv.rates
        ]
    )

else:
    print(
        "  gamma categories: 1"
    )


print("\nComposition components and resulting Q matrices:")
print("Q row/column order: A C G T")


q_matrices = {}


def dna_q_matrix(comp_values, r_values):
    """
    Build a normalized DNA Q matrix from p4's ACGT frequencies and
    AC, AG, AT, CG, CT, GT exchangeabilities.
    """
    matrix = [
        [0.0 for column in range(4)]
        for row in range(4)
    ]

    rate_pairs = [
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
        (2, 3)
    ]

    for rate, pair in zip(r_values, rate_pairs):
        left, right = pair

        matrix[left][right] = float(rate) * float(comp_values[right])
        matrix[right][left] = float(rate) * float(comp_values[left])

    for row in range(4):
        matrix[row][row] = -sum(
            matrix[row][column]
            for column in range(4)
            if column != row
        )

    expected_rate = -sum(
        float(comp_values[row]) * matrix[row][row]
        for row in range(4)
    )

    if expected_rate <= 0.0:
        raise RuntimeError(
            "Cannot normalize Q matrix with non-positive expected rate."
        )

    return [
        [
            value / expected_rate
            for value in row
        ]
        for row in matrix
    ]


for comp in part.comps:
    q_matrix = dna_q_matrix(
        comp.val,
        r_matrix.val
    )

    q_matrices[int(comp.num)] = q_matrix

    print(
        "\nComponent %d: free=%d"
        % (
            display_comp_num(comp.num),
            int(comp.free)
        )
    )

    print(
        "Frequencies: %s"
        % [
            float(value)
            for value in comp.val
        ]
    )

    print("Resulting Q:")
    print("             A            C            G            T")

    for state, row in zip("ACGT", q_matrix):
        print(
            "%s  %12.8f %12.8f %12.8f %12.8f"
            % (
                (state,) + tuple(row)
            )
        )


# ============================================================
# Print branch assignments
# ============================================================

print("\nBranch composition assignments:")

for node in t.iterNodesNoRoot():
    node_number = int(node.nodeNum)

    if node_number in selected_branch_records:
        selected_by = format_selection_sources(
            selected_branch_records[node_number]["sources"]
        )
    else:
        selected_by = "background"

    print(
        "%s : comp %d, length=%.10g, descendants=%s, selected_by=%s"
        % (
            branch_label(node),
            display_comp_num(node.parts[0].compNum),
            float(node.br.len),
            ",".join(descendant_taxa(node)),
            selected_by
        )
    )


# ============================================================
# Save optimized trees
# ============================================================

t.writeNewick(TREE_OUTPUT_FILE)


optimized_branch_lengths = {}

for node in t.iterNodesNoRoot():
    optimized_branch_lengths[int(node.nodeNum)] = float(node.br.len)

    node.br.len = float(
        display_comp_num(node.parts[0].compNum)
    )


try:
    t.writeNewick(COMPONENT_TREE_OUTPUT_FILE)

finally:
    for node in t.iterNodesNoRoot():
        node.br.len = optimized_branch_lengths[
            int(node.nodeNum)
        ]


# ============================================================
# Construct JSON output
# ============================================================

result = {
    "logL": float(t.logLike),

    "input": {
        "alignment": ALIGNMENT_FILE,
        "tree": TREE_FILE
    },

    "output_prefix": args.pre,

    "requested_comp_2_paths": requested_specs,

    "requested_comp_2_branches": requested_specs,

    "requested_comp_2_clades": requested_clades,

    "number_of_unique_comp_2_branches": len(
        selected_branch_records
    ),

    "composition_optimization": {
        "equal_frequencies_fixed": bool(
            args.equal_comp_freq
        ),

        "frequencies_optimized": bool(
            not args.equal_comp_freq
        ),

        "optimizer": args.optimizer,

        "refined_with_newtAndBrentPowell": bool(
            refine_likelihood
        ),

        "p4_number_of_free_model_parameters": int(
            t.model.nFreePrams
        )
    },

    "resolved_comp_2_paths": [],

    "resolved_comp_2_clades": [],

    "substitution_model": {
        "name": model_name,

        "requested_name": model_definition["requested_name"],

        "base_model": base_model,

        "iqtree_rate_code": model_definition["rate_code"],

        "rate_mode": model_definition["rate_mode"],

        "equal_frequencies_in_iqtree": bool(
            model_definition["equal_frequencies_in_iqtree"]
        ),

        "ndch_frequencies_optimized": bool(
            not args.equal_comp_freq
        ),

        "equal_comp_frequencies_fixed": bool(
            args.equal_comp_freq
        ),

        "invariant_sites": {
            "proportion": float(part.pInvar.val),
            "free": int(part.pInvar.free)
        },

        "r_matrix": int(r_matrix.num),

        "spec": r_matrix.spec,

        "free": int(r_matrix.free),

        "parameter_names": (
            ["kappa"]
            if r_matrix.spec == "2p"
            else ["AC", "AG", "AT", "CG", "CT", "GT"]
        ),

        "values": [
            float(value)
            for value in r_matrix.val
        ],

        "gamma": (
            {
                "n_categories": int(part.nGammaCat),

                "alpha": float(gdasrv.val[0]),

                "free": int(gdasrv.free),

                "category_frequencies": [
                    float(value)
                    for value in gdasrv.freqs
                ],

                "category_rates": [
                    float(value)
                    for value in gdasrv.rates
                ]
            }
            if gdasrv is not None
            else None
        )
    },

    "q_state_order": list("ACGT"),

    "compositions": {},

    "q_matrices": {},

    "branches": {}
}


for resolved_path in resolved_paths:
    result["resolved_comp_2_paths"].append(
        {
            "specification": resolved_path["specification"],

            "left_group": resolved_path["left_text"],

            "right_group": resolved_path["right_text"],

            "left_taxa": resolved_path["left_taxa"],

            "right_taxa": resolved_path["right_taxa"],

            "left_lca": node_label(
                resolved_path["left_lca"]
            ),

            "right_lca": node_label(
                resolved_path["right_lca"]
            ),

            "number_of_branches": len(
                resolved_path["branches"]
            ),

            "branches": [
                {
                    "parent": node_label(node.parent),

                    "child": node_label(node),

                    "node_number": int(node.nodeNum),

                    "descendant_taxa": descendant_taxa(node)
                }
                for node in resolved_path["branches"]
            ]
        }
    )


for resolved_clade in resolved_clades:
    result["resolved_comp_2_clades"].append(
        {
            "specification": resolved_clade["specification"],

            "taxa": resolved_clade["taxa"],

            "lca": node_label(
                resolved_clade["lca"]
            ),

            "lca_descendant_taxa": descendant_taxa(
                resolved_clade["lca"]
            ),

            "number_of_branches": len(
                resolved_clade["branches"]
            ),

            "branches": [
                {
                    "parent": node_label(node.parent),

                    "child": node_label(node),

                    "node_number": int(node.nodeNum),

                    "descendant_taxa": descendant_taxa(node)
                }
                for node in resolved_clade["branches"]
            ]
        }
    )


for comp in part.comps:
    internal_number = int(comp.num)
    output_number = display_comp_num(internal_number)

    start_values = starting_compositions[internal_number]
    end_values = optimized_compositions[internal_number]

    result["compositions"][str(output_number)] = {
        "free": int(comp.free),

        "starting_values": [
            float(value)
            for value in start_values
        ],

        "optimized_values": [
            float(value)
            for value in end_values
        ],

        "maximum_absolute_change": float(
            maximum_absolute_change(
                start_values,
                end_values
            )
        )
    }

    result["q_matrices"][str(output_number)] = [
        [
            float(value)
            for value in row
        ]
        for row in q_matrices[internal_number]
    ]


for node in t.iterNodesNoRoot():
    node_number = int(node.nodeNum)

    if node_number in selected_branch_records:
        selection_sources = list(
            selected_branch_records[node_number]["sources"]
        )
    else:
        selection_sources = []

    result["branches"][node_label(node)] = {
        "node_number": node_number,

        "parent": node_label(node.parent),

        "length": float(node.br.len),

        "comp": display_comp_num(
            node.parts[0].compNum
        ),

        "descendant_taxa": descendant_taxa(node),

        "selection_sources": selection_sources
    }


with open(RESULT_FILE, "w") as output_file:
    json.dump(
        result,
        output_file,
        indent=2
    )


# ============================================================
# Human-readable report
# ============================================================

def read_text_file(file_name):
    """
    Read and strip a text file.
    """
    with open(file_name) as input_file:
        return input_file.read().strip()


with open(REPORT_FILE, "w") as report:
    report.write("P4 NDCH ANALYSIS REPORT\n")
    report.write("======================\n\n")

    report.write(
        "Input alignment: %s\n"
        % ALIGNMENT_FILE
    )

    report.write(
        "Input tree:      %s\n"
        % TREE_FILE
    )

    report.write(
        "Requested model: %s\n"
        % model_name
    )

    report.write(
        "Log likelihood:  %.10f\n"
        % float(t.logLike)
    )

    report.write(
        "Optimizer:       %s\n"
        % args.optimizer
    )

    report.write(
        "Refinement:      %s\n"
        % (
            "newtAndBrentPowell"
            if refine_likelihood
            else "none"
        )
    )

    report.write(
        "Free parameters: %d\n"
        % int(t.model.nFreePrams)
    )

    report.write(
        "Composition frequencies: %s\n"
        % (
            "fixed equal"
            if args.equal_comp_freq
            else "independently optimized"
        )
    )

    report.write(
        "Requested paths: %s\n"
        % (
            "; ".join(requested_specs)
            if requested_specs
            else "none"
        )
    )

    report.write(
        "Requested clades: %s\n"
        % (
            "; ".join(requested_clades)
            if requested_clades
            else "none"
        )
    )

    report.write(
        "Unique component 2 branches: %d\n\n"
        % len(selected_branch_records)
    )

    report.write("COMPOSITION OPTIMIZATION\n")
    report.write("------------------------\n\n")

    for comp in part.comps:
        internal_number = int(comp.num)
        output_number = display_comp_num(internal_number)

        start_values = starting_compositions[internal_number]
        end_values = optimized_compositions[internal_number]

        report.write(
            "Component %d\n"
            % output_number
        )

        report.write(
            "  free: %d\n"
            % int(comp.free)
        )

        report.write(
            "  starting A C G T: %s\n"
            % " ".join(
                "%.12g" % value
                for value in start_values
            )
        )

        report.write(
            "  optimized A C G T: %s\n"
            % " ".join(
                "%.12g" % value
                for value in end_values
            )
        )

        report.write(
            "  maximum absolute change: %.12g\n\n"
            % maximum_absolute_change(
                start_values,
                end_values
            )
        )

    report.write("RESOLVED COMPONENT 2 PATHS\n")
    report.write("--------------------------\n\n")

    if resolved_paths:
        for resolved_path in resolved_paths:
            report.write(
                "Path: %s\n"
                % resolved_path["specification"]
            )

            report.write(
                "  Left: LCA(%s) = %s\n"
                % (
                    resolved_path["left_text"],
                    node_label(resolved_path["left_lca"])
                )
            )

            report.write(
                "  Right: LCA(%s) = %s\n"
                % (
                    resolved_path["right_text"],
                    node_label(resolved_path["right_lca"])
                )
            )

            for index, node in enumerate(
                resolved_path["branches"],
                start=1
            ):
                report.write(
                    "    %d. %s; descendants=%s\n"
                    % (
                        index,
                        branch_label(node),
                        ",".join(descendant_taxa(node))
                    )
                )

            report.write("\n")

    else:
        report.write("none\n\n")

    report.write("RESOLVED COMPONENT 2 CLADES\n")
    report.write("---------------------------\n\n")

    if resolved_clades:
        for resolved_clade in resolved_clades:
            report.write(
                "Clade: %s\n"
                % resolved_clade["specification"]
            )

            report.write(
                "  LCA: %s\n"
                % node_label(resolved_clade["lca"])
            )

            report.write(
                "  LCA descendants: %s\n"
                % ",".join(
                    descendant_taxa(resolved_clade["lca"])
                )
            )

            for index, node in enumerate(
                resolved_clade["branches"],
                start=1
            ):
                report.write(
                    "    %d. %s; descendants=%s\n"
                    % (
                        index,
                        branch_label(node),
                        ",".join(descendant_taxa(node))
                    )
                )

            report.write("\n")

    else:
        report.write("none\n\n")

    report.write("SUBSTITUTION PROCESS\n")
    report.write("--------------------\n\n")

    report.write(
        "Rate matrix specification: %s\n"
        % r_matrix.spec
    )

    report.write(
        "Rate matrix free: %d\n"
        % int(r_matrix.free)
    )

    report.write(
        "Rate values: %s\n"
        % " ".join(
            "%.12g" % float(value)
            for value in r_matrix.val
        )
    )

    report.write(
        "Invariant proportion: %.12g; free=%d\n"
        % (
            float(part.pInvar.val),
            int(part.pInvar.free)
        )
    )

    if gdasrv is not None:
        report.write(
            "Gamma categories: %d\n"
            % int(part.nGammaCat)
        )

        report.write(
            "Gamma alpha: %.12g; free=%d\n"
            % (
                float(gdasrv.val[0]),
                int(gdasrv.free)
            )
        )

        report.write(
            "Gamma rates: %s\n"
            % " ".join(
                "%.12g" % float(value)
                for value in gdasrv.rates
            )
        )

    else:
        report.write(
            "Gamma categories: 1\n"
        )

    report.write("\nQ MATRICES\n")
    report.write("----------\n")

    for comp in part.comps:
        internal_number = int(comp.num)

        report.write(
            "\nComponent %d; rows/columns A C G T\n"
            % display_comp_num(internal_number)
        )

        for state, row in zip(
            "ACGT",
            q_matrices[internal_number]
        ):
            report.write(
                "  %s  %s\n"
                % (
                    state,
                    " ".join(
                        "%12.8f" % float(value)
                        for value in row
                    )
                )
            )

    report.write("\nMAXIMUM-LIKELIHOOD TREE\n")
    report.write("-----------------------\n\n")

    report.write(
        "Optimized branch lengths:\n%s\n\n"
        % read_text_file(TREE_OUTPUT_FILE)
    )

    report.write(
        "Components encoded as branch lengths:\n%s\n\n"
        % read_text_file(COMPONENT_TREE_OUTPUT_FILE)
    )

    report.write("BRANCH COMPONENT ASSIGNMENTS\n")
    report.write("----------------------------\n\n")

    report.write(
        "parent -> child  component  length  "
        "descendant_taxa  selected_by\n"
    )

    for node in t.iterNodesNoRoot():
        node_number = int(node.nodeNum)

        if node_number in selected_branch_records:
            selected_by = format_selection_sources(
                selected_branch_records[node_number]["sources"]
            )
        else:
            selected_by = "background"

        report.write(
            "%s  %d  %.12g  %s  %s\n"
            % (
                branch_label(node),
                display_comp_num(
                    node.parts[0].compNum
                ),
                float(node.br.len),
                ",".join(descendant_taxa(node)),
                selected_by
            )
        )


# ============================================================
# Final summary
# ============================================================

print("\nFinal optimized composition vectors:")

for comp in part.comps:
    print(
        "  Component %d: free=%d, ACGT=%s"
        % (
            display_comp_num(comp.num),
            int(comp.free),
            [
                float(value)
                for value in comp.val
            ]
        )
    )


print("\nAnalysis results written to:")

print(
    "  Maximum-likelihood tree: %s"
    % TREE_OUTPUT_FILE
)

print(
    "  Component-labeled tree: %s"
    % COMPONENT_TREE_OUTPUT_FILE
)

print(
    "  JSON results:            %s"
    % RESULT_FILE
)

print(
    "  NDCH report:             %s"
    % REPORT_FILE
)
