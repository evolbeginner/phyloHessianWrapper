#! /usr/bin/env Rscript


##################################
suppressPackageStartupMessages(
    suppressWarnings({
    library(phytools)
    library(phangorn)
    })
)


##################################
args <- commandArgs(trailingOnly = T)

treefile <- args[1]
phy <- read.tree(treefile)

branch_out_infile <- args[2]
df <- read.table(branch_out_infile, header = F)
colnames(df) <- c('branch', 'bl')
df <- cbind(df, "tree_order" = seq_len(nrow(df)))

nb.tip <- length(phy$tip.label)
nb.node <- phy$Nnode


##################################
branches <- vector()
j <- 0
all_children <- Children(phy, (1 + nb.tip):(nb.node + nb.tip) )

for(anc in (nb.node+nb.tip):(1+nb.tip)){
    #if(anc == nb.tip+1){next}
    children <- all_children[[anc-nb.tip]]

    for(i in 1:length(children)){
        j <- j + 1L
        child <- children[i]
        tips <- Descendants(phy, child, type='tip')
        tip_names <- sort(as.vector(sapply(tips, function(tip){phy$tip.label[tip]})))
        tip_name <- paste(tip_names, collapse='-')
        comp_tip_names <- sort( phy$tip.label[ !(phy$tip.label %in% tip_names) ] )
        comp_tip_name <- paste(comp_tip_names, collapse='-')
        branch <- paste(sort(c(tip_name, comp_tip_name)), collapse=',')
        branches <- c(branches, branch)
    }
}

#df$ape_order <- sapply(1:nrow(df), function(i){ which(branches == df[i,]$branch) })
df$ape_order <- sapply(
    1:nrow(df), function(i){
        alt_branch <- paste(rev(strsplit(df[i,]$branch, ",")[[1]]), collapse = ",")
        which(branches == df[i,]$branch | branches == alt_branch)
        #which(branches == df[i,]$branch)
    }
)

write.table(df[order(df[,1]),], sep="\t", quote=F, row.names=T);

cat(paste(df[order(df$ape_order),]$bl, collapse=','))
cat("\n")


