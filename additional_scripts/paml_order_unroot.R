#!/usr/bin/env Rscript

# paml_order_unroot.R
# Convert ONE rooted binary Newick tree to PAML-style unrooted ordering,
# output topology only (no branch lengths).
#
# Usage:
#   Rscript paml_order_unroot.R input_file output.nwk
#
# Optional:
#   SPLIT_CHILD=1 or 2
#   - If root has two internal children, choose which root child to split.
#
# Notes:
# - Supports species.trees-like files:
#     line 1: metadata (e.g., "5 1")
#     line 2: rooted tree (possibly with trailing tags like >8<12;)
# - Auto-detects tree in line 1 or 2 (then scans all lines if needed).

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript paml_order_unroot.R input_file output.nwk")
}

infile  <- args[1]
outfile <- args[2]

raw_lines <- readLines(infile, warn = FALSE)
lines <- trimws(raw_lines)
lines <- lines[nzchar(lines)]

if (length(lines) == 0) stop("Input file has no non-empty lines.")

clean_newick_line <- function(x) {
  x <- trimws(x)
  # Remove trailing PAML-like annotations such as >8<12 before ';'
  # Example: (t5,((t4,t1),(t3,t2)))>8<12; -> (t5,((t4,t1),(t3,t2)));
  x <- sub(">[^;]*", "", x)
  if (!grepl(";", x, fixed = TRUE)) x <- paste0(x, ";")
  x
}

looks_like_newick <- function(x) {
  # Metadata lines like "5 1" should fail this quickly
  grepl("\\(", x) && grepl("\\)", x)
}

try_parse <- function(x) {
  x2 <- clean_newick_line(x)
  tryCatch(read.tree(text = x2), error = function(e) NULL)
}

is_rooted_binary <- function(tr) {
  if (is.null(tr)) return(FALSE)
  if (length(tr$tip.label) < 2) return(FALSE)
  isTRUE(is.rooted(tr)) && isTRUE(is.binary(tr))
}

# Detection order: line 1, line 2 first, then the rest
order_idx <- unique(c(1, 2, seq_along(lines)))
order_idx <- order_idx[order_idx >= 1 & order_idx <= length(lines)]

tr <- NULL
chosen_line <- NA_integer_

for (i in order_idx) {
  ln <- lines[i]
  if (!looks_like_newick(ln)) next
  tmp <- try_parse(ln)
  if (is_rooted_binary(tmp)) {
    tr <- tmp
    chosen_line <- i
    break
  }
}

if (is.null(tr)) {
  stop("Could not find a rooted binary Newick tree in the input file.")
}

message(sprintf("Detected rooted binary tree on non-empty line %d", chosen_line))

Ntip <- length(tr$tip.label)
root <- Ntip + 1
edge <- tr$edge

children_of <- function(node) edge[edge[, 1] == node, 2]

emit_subtree_topology <- function(node) {
  if (node <= Ntip) return(tr$tip.label[node])
  kids <- children_of(node)
  parts <- vapply(kids, emit_subtree_topology, character(1))
  paste0("(", paste(parts, collapse = ", "), ")")
}

root_kids <- children_of(root)
if (length(root_kids) != 2) stop("Root must have exactly 2 children.")

# Choose which root child to split
split_env <- Sys.getenv("SPLIT_CHILD", unset = "")
if (split_env %in% c("1", "2")) {
  split_child <- root_kids[as.integer(split_env)]
} else {
  internal <- root_kids[root_kids > Ntip]
  if (length(internal) == 1) {
    split_child <- internal[1]   # split internal side if other side is a tip
  } else {
    split_child <- root_kids[1]  # fallback
  }
}
other_child <- setdiff(root_kids, split_child)

if (split_child <= Ntip) {
  stop("Chosen split_child is a tip. Use SPLIT_CHILD=1 or 2 to choose the other side.")
}

split_kids <- children_of(split_child)
if (length(split_kids) != 2) stop("split_child must have exactly 2 children.")

k1 <- split_kids[1]
k2 <- split_kids[2]

# PAML-style unrooted top-level order: (child_of_split_1, child_of_split_2, other_side)
out_newick <- paste0(
  "(",
  emit_subtree_topology(k1), ", ",
  emit_subtree_topology(k2), ", ",
  emit_subtree_topology(other_child),
  ");"
)

writeLines(out_newick, outfile)
cat(out_newick, "\n")
