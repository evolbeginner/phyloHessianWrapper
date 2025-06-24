#! /usr/bin/env Rscript


##################################################
suppressPackageStartupMessages(
    suppressWarnings({
    library(phangorn)
    library(phytools)
    library(getopt)
    })
)


##################################################
spec <- matrix(c(
	'tree', 't', 2, 'character'
    ), byrow=T, ncol=4
)
opts <- getopt(spec)

if(!is.null(opts$tree)){
    t <- read.tree(opts$tree)
}

nb.tip <- length(t$tip.label)
nb.node <- t$Nnode


##################################################
#t$edge.length[t$edge.length < 0.1] <- 0.1

all_children <- Children(t, (1 + nb.tip):(nb.node + nb.tip) )

output <- vector()
for(i in length(all_children):1){
    #sapply(all_children[[i]], function(x){ cat(i+nb.tip, x, t$edge.length[t$edge[,2] == x], "\n") })
    #cat(paste0(paste0(sapply(all_children[[i]], function(x){ t$edge.length[t$edge[,2] == x] }), collapse=','), ','))
    output <- c(output, paste0(sapply(all_children[[i]], function(x){ t$edge.length[t$edge[,2] == x] }), collapse=','))
}
cat(paste(output, collapse=','), "\n")


