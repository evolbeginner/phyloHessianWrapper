#! /usr/bin/env Rscript

######################################
# Update
# 2025-02-06
#   make seq file (-s) and treefile phy$tip.label the same order

set.seed(0)

######################################
suppressPackageStartupMessages(
suppressWarnings({
    library(getopt)
    #library(parallel)
    #library(matrixStats)
    library(seqinr)
    #library(phytools)
    library(phangorn)
})
)

######################################
DIR <- commandArgs(trailingOnly = FALSE)
DIR <- paste0(dirname(sub("--file=", "", DIR[grepl("--file=", DIR)])), '/')
DIR <- paste0(normalizePath(DIR), '/') # absolute path

E <- ape::matexpo

AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-")
AA_indices <- c(1:20, 999)
AA_list <- setNames(AA_indices, AAs)

format <- 'fasta'

######################################
# Helper Functions
read_traits <- function(d){
    v <- list()
    name <- unlist(d$V1)
    for(i in 2:ncol(d)){
        v[[i-1]] <- unlist(d[,i])
        names(v[[i-1]]) <- name
    }
    return(v)
}

rename_trait <- function(x){
    if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label))
            x <- x[phy$tip.label]
        else stop("the names of 'x' and the tip labels of the tree do not match")
    }
    if (!is.factor(x)) {
        x <- factor(x)
        x <- as.integer(x)
    }
    x
}

get_v_freq <- function(v){
    u <- lapply(v, function(x){paste(x, collapse='!')})
    t <- table(unlist(u))
    t_names_split <- lapply(names(t), function(x){
        unlist(strsplit(x,"!"), use.name=F)
    })
    l <- vector("list", length(rownames(t)))
    for(i in 1:length(rownames(t))){
        l[[i]] <- list(as.integer(t_names_split[[i]]), t[[i]])
    }
    l
}

replace_not_in_AAs <- function(x) {
    x[is.na(x)] <- 999
    return(x)
}

deduplicate_pattern <- function(v){
    v2 <- list()
    site <- numeric(length(v))
    duplicate_count <- setNames(integer(0), character(0))

    for (i in seq_along(v)) {
        element <- v[[i]][1]
        element_str <- paste(element[[1]], collapse = ",")
        if (!element_str %in% names(duplicate_count)) {
            duplicate_count[element_str] <- 1
        } else {
            duplicate_count[element_str] <- duplicate_count[element_str] + 1
        }
        site[i] <- which(names(duplicate_count) == element_str)
    }
    
    for (name in names(duplicate_count)) {
        element_vector <- strsplit(name, ",")[[1]]
        v2[[length(v2) + 1]] <- list(
            as.character(element_vector), 
            unname(duplicate_count[name])
        )
    }
    return(list(v2=v2, site=site))
}

find_cherry_nodes <- function(phy) {
    cherry_nodes <- vector()
    for (node in ((length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode))) {
        children <- all_children[node - length(phy$tip.label)]
        is_tip_children <- sapply(children, function(node){
            node <= length(phy$tip.label)
        })
        if (all(is_tip_children)) cherry_nodes <- append(cherry_nodes, node)
    }
    return(cherry_nodes)
}

output_julia <- function(julia_outdir, all_children, v){
    if(!is.null(julia_outdir)){
        # Output basics
        outfile <- file.path(julia_outdir,"basics")
        out_fh <- file(outfile, "w")
        cat(paste("nb.node", nb.node, sep="\t"), "\n", file=out_fh, sep="")
        cat(paste("nb.tip", nb.tip, sep="\t"), "\n", file=out_fh, sep="")
        close(out_fh)

        # Output all_children
        outfile <- file.path(julia_outdir,"all_children")
        out_fh <- file(outfile, "w")
        sapply(1:length(all_children), function(i){
            cat(paste(unlist(c(i+length(phy$tip.label),all_children[[i]])), 
                collapse="\t"), "\n", file=out_fh, sep='')
        })
        close(out_fh)

        # Output descendants
        outfile <- file.path(julia_outdir,"descendants")
        out_fh <- file(outfile, "w")
        sapply((length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode), 
            function(i){
                cat(paste(unlist(c(i, Descendants(phy,i,type="tips"))),
                    collapse="\t"), "\n", file=out_fh, sep="")
            })
        close(out_fh)

        # Output all descendants
        outfile <- file.path(julia_outdir,"all")
        out_fh <- file(outfile, "w")
        sapply((length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode), 
            function(i){
                cat(paste(unlist(c(i, Descendants(phy,i,type="all"))),
                    collapse="\t"), "\n", file=out_fh, sep="")
            })
        close(out_fh)

        # Output pattern
        outfile <- file.path(julia_outdir,"pattern")
        out_fh <- file(outfile, "w")
        sapply(1:length(v), function(i){
            cat(paste(unlist(c(v[[i]][[1]],v[[i]][[2]])),
                collapse="\t"), "\n", file=out_fh, sep="")
        })
        close(out_fh)

        # Output cherry nodes
        outfile <- file.path(julia_outdir,"cherry")
        out_fh <- file(outfile, "w")
        cat(paste(cherry_nodes,collapse="\n"), "\n", file=out_fh, sep="")
        close(out_fh)

        # Output site2pattern
        outfile <- file.path(julia_outdir,"site2pattern")
        out_fh <- file(outfile, 'w')
        output <- sapply(1:length(site), function(i){
            paste(i,site[i],sep="\t")
        })
        write(output, file=out_fh)
        close(out_fh)
    }
    else{
        return("no need to create julia_outdir")
    }
}

######################################
# Command line parsing
spec <- matrix(c(
    'tree', 't', 2, 'character',
    'format', 'F', 2, 'character',
    'traits', 's', 2, 'character',
    'cpu', 'n', 2, 'integer',
    'type', '', 2, 'character',
    'julia_outdir', 'j', 2, 'character',
    'force', 'f', 0, 'logical'
    ), ncol=4, byrow=T
)

opts <- getopt(spec)

######################################
# Main execution
if(!is.null(opts$tree)){
    phy <- read.tree(opts$tree)
    phy <- unroot(phy)
}

if(!is.null(opts$julia_outdir)){
    julia_outdir <- opts$julia_outdir
    if(!dir.exists(julia_outdir)){
        dir.create(julia_outdir, recursive=T)
    } else if(opts$force){
        unlink(julia_outdir, recursive=T)
        dir.create(julia_outdir, recursive=T)
    } else{
        stop(paste("julia_outdir", julia_outdir, "already exists! Exiting ......"))
    }
}

# Process sequence data
if(!is.null(opts$traits)){
    s <- read.fasta(opts$traits)
    if(var(getLength(s)) != 0){
        stop("seqfile: seq of diff lengths")
    }
    seq <- getSequence(s)

    v <- list()
    for(i in 1:getLength(s)[1]){
        v[[i]] <- sapply(s, function(x) x[i])
        v[[i]] <- v[[i]][match(phy$tip.label, names(v[[i]]))]
    }

    nl <- nlevels(factor(unlist(v, use.names=F)))

    if(opts$type == "DNA"){
        v <- lapply(v, rename_trait)
        cat("Alignment length:", length(v), "\n")
        v <- get_v_freq(v)
        cat("No. of pattern:", length(v), "\n")
    } else if(opts$type == "AA"){
        v <- lapply(v, function(x) unname(AA_list[toupper(x)]))
        for (i in seq_along(v)) {
            v[[i]] <- replace_not_in_AAs(v[[i]])
        }
        v <- lapply(v, function(x) list(x,1))
        dedup_rv <- deduplicate_pattern(v)
        v <- dedup_rv$v2
        site <- dedup_rv$site
    } else {
        stop("unknown seq type! Exiting ......")
    }
}

# Calculate tree information
nb.tip <- length(phy$tip.label)
nb.node <- phy$Nnode
all_children <- Children(phy, (1 + nb.tip):(nb.node + nb.tip))
cherry_nodes <- find_cherry_nodes(phy)

# Output to Julia
output_julia(julia_outdir, all_children, v)
