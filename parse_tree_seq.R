#! /bin/env Rscript


######################################
# Update
# 2025-02-06
#   make seq file (-s) and treefile phy$tip.label the same order


######################################
# multiple Qs
# EX2+F now working
# +G+I works for EX2+G+I

# note EX2.dat0 should be used
# also is_inv_site prev wrong; weight not integrated

# working on PMSF


######################################
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
DIR <- paste0( dirname(sub("--file=", "", DIR[grepl("--file=", DIR)])), '/' )
DIR <- paste0(normalizePath(DIR), '/') # absolute path

E <- ape::matexpo

AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-")
AA_indices <- c(1:20, 999)
AA_list <- setNames(AA_indices, AAs)

DNA_list <- setNames(c(1:4, 999), c("A", "C", "G", "T", "-"))

format <- 'fasta'


######################################
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
		else stop("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
	}

	x <- toupper(x)
	x <- unname(DNA_list[x])
	x[is.na(x)] <- 999
	x
}


get_v_freq <- function(v){
	u <- lapply(v, function(x){paste(x, collapse='!')})
	t <- table(unlist(u))
	t_names_split <- lapply(names(t), function(x){unlist(strsplit(x,"!"), use.name=F)})
	l <- vector("list", length(rownames(t)))
	for(i in 1:length(rownames(t))){
		l[[i]] <- list(as.integer(t_names_split[[i]]), t[[i]])
	}
	l
}


deduplicate_pattern <- function(v){
    v2 <- list()
    site <- numeric(length(v))

    # Create a named vector to keep track of duplicate counts
    duplicate_count <- setNames(integer(0), character(0))

    for (i in seq_along(v)) {
        element <- v[[i]][1]
        element_str <- paste(element[[1]], collapse = ",")  # Convert the vector to a string
        # Update the count of the current element
        if (!element_str %in% names(duplicate_count)) {
            duplicate_count[element_str] <- 1
        } else {
            duplicate_count[element_str] <- duplicate_count[element_str] + 1
        }
        site[i] <- which(names(duplicate_count) == element_str)
    }
    # Create v2 with unique elements and their counts
    for (name in names(duplicate_count)) {
        element_vector <- strsplit(name, ",")[[1]]  # Convert the string back to a vector
        v2[[length(v2) + 1]] <- list(as.character(element_vector), unname(duplicate_count[name]))
    }
    return(list(v2=v2, site=site))
}


replace_not_in_AAs <- function(x) {
  x[is.na(x)] <- 999
  return(x)
}


######################################
output_julia <- function(julia_outdir, all_children, v){
    if(!is.null(julia_outdir)){
        outfile <- file.path(julia_outdir,"basics"); out_fh <- file(outfile, "w")
        cat(paste("nb.node", nb.node, sep="\t"), "\n", file=out_fh, sep="")
        cat(paste("nb.tip", nb.tip, sep="\t"), "\n", file=out_fh, sep="")
        close(out_fh)

        outfile <- file.path(julia_outdir,"all_children"); out_fh <- file(outfile, "w")
        output_list <- sapply(1:length(all_children), function(i){cat(paste(unlist(c(i+length(phy$tip.label),all_children[[i]])), collapse="\t"), "\n", file=out_fh, sep='')})
        close(out_fh)

        outfile <- file.path(julia_outdir,"descendants"); out_fh <- file(outfile, "w")
        output_list <- sapply((length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode), function(i){cat(paste(unlist(c(i, Descendants(phy,i,type="tips"))),collapse="\t"), "\n", file=out_fh, sep="")})
        close(out_fh)

        outfile <- file.path(julia_outdir,"all"); out_fh <- file(outfile, "w")
        output_list <- sapply((length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode), function(i){cat(paste(unlist(c(i, Descendants(phy,i,type="all"))),collapse="\t"), "\n", file=out_fh, sep="")})
        close(out_fh)

        outfile <- file.path(julia_outdir,"pattern"); out_fh <- file(outfile, "w")
        output_list <- sapply(1:length(v), function(i){cat(paste(unlist(c(v[[i]][[1]],v[[i]][[2]])),collapse="\t"), "\n", file=out_fh, sep="")})
        close(out_fh)

        outfile <- file.path(julia_outdir,"cherry"); out_fh <- file(outfile, "w")
        cat(paste(cherry_nodes,collapse="\n"), "\n", file=out_fh, sep="")
        close(out_fh)

        outfile <- file.path(julia_outdir,"site2pattern"); out_fh <- file(outfile, 'w')
        output <- sapply( 1:length(site), function(i){paste(i,site[i],sep="\t")} )
        write(output, file=out_fh)
        close(out_fh)
    }
    else{
        return(paste("no need to create julia_outdir"))
    }
}


######################################
find_cherry_nodes <- function(phy) {
    cherry_nodes <- vector()
    for (node in ((length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode))) {
        children <- all_children[node - length(phy$tip.label)]
        is_tip_children <- sapply(children, function(node){node <= length(phy$tip.label)})
        if (all(is_tip_children)) cherry_nodes <- append(cherry_nodes, node)
    }
    return(cherry_nodes)
}


######################################
do_phylo_log_lk <- function(param, is_inv=F, is_inv_site=F) {
	liks_ori <- liks
	bl <- param
    rs <- c(1)
	#rs <- c(0.6086, 2.9268)
    #rs <- c(0.0207, 1.9793)
    #print(Qs[[1]]); q()

    is_inv <- F
	if(is_inv){
		inv <- tail(param, n=1)
        # for a given inv
        inv <- 38.71/100
		#rs <- rs/(1-inv)
        rs <- c(0, rs)
	}

	log_lk <- 0
	lk <- numeric(length(rs))
    counter <- 0

    # rate_fixed for Q (Qrs)
    Qrs <- rep(1, length(Qs))
    #Qrs <- c(.4448, 1.1213)
    #Qrs <- c(.8645, 2.1115) # 4-taxon 2.aln
    #Qrs <- c(.7452, 1.8201)
    #Qrs <- c(0.6517, 1.5917) #EX2
    #Qrs <- c(0.557, 1.360)
    #Qrs <- c(0.5044, 0.9880, 1.7913) #EX3

    prop <- rep(1/length(rs), length(rs))
    freq <- rep(1/length(Qrs), length(Qrs))
    #freq <- c(0.2, 0.8) # freq is the "raw weight" of each Q, while "weight" is the final weight
    #freq <- c(0.682, 0.318)
    #freq <- c(0.448, 0.552)
    #freq <- c(0.763, 0.237)
    #freq <- c(0.1793, 0.8207)
    #freq <- c(0.6134, 0.3866)
    #freq <- c(0.891, .109) # 4-taxon 2.aln
    weight <- outer(prop, freq, "*"); weight <- c(t(weight))
    #Qrs <- sapply( Qrs, function(x){x/sum(Qrs*freq)} )
    #print(freq * sapply(1:length(Qs), function(i){sum( Pis[[i]] * diag(Qs[[i]]) * Qrs[i] )})); q()
    #Qs <- lapply( Qs, function(x){x/sum(Qrs*freq)} )
    #print(freq * sapply(1:length(Qs), function(i){sum( Pis[[i]] * diag(Qs[[i]]) * Qrs[i] )})); q()
    #print(lapply(Qs, diag)); q()
    #print( sum( freq[1]*Qrs[1]*diag(Qs[[1]]) + freq[2]*Qrs[2]*diag(Qs[[2]]) ) )
    #print(freq * Qrs); q()

	for(r_i in 1:length(rs)){
        for(q_i in 1:length(Qs)){
            counter <- counter + 1
            # define the Q and Pi
            Q <- Qs[[q_i]]; Pi <- Pis[[q_i]]

            # rate_fixed for Q (Qrs)
            r <- Qrs[q_i] * rs[r_i]

            j <- 0
            liks <- liks_ori
            comp <- numeric(nb.tip + nb.node)

            for (anc in (nb.node + nb.tip):(1 + nb.tip)) {
                children <- all_children[[anc-nb.tip]]
                m <- matrix(0, nl, length(children))
                if(as.character(anc) %in% names(cherry_liks)){
                    l <- cherry_liks[[as.character(anc)]]
                    site_pattern <- paste(apply(liks[children,], 1, paste, collapse = ""), collapse = "")
                    m_prod <- l[[site_pattern]]
                    j <- j + length(children)
                }else{
                    for(i in 1:length(children)){
                        j <- j + 1L
                        #cat(paste(c(anc,j),collapse="|"), children[i], bl[j], paste(unlist(Descendants(phy,j,type='tips')),collapse='-'), "\n")
                        m[,i] <- E(Q * (bl[j])*r) %*% liks[children[i], ]
                        #m[,i] <- U %*% diag(exp(Lambda * bl[j]*r )) %*% U_inv %*% liks[children[i], ]
                        #print(c(anc, children[i]))
                    }
                    m_prod <- rowProds(m)
                }
                comp[anc] <- sum(m_prod)
                #if(anc == (1 + nb.tip)){ comp[anc] <- rep(1/nl,nl) %*% m_prod }
                if(anc == (1 + nb.tip)){ comp[anc] <- Pi %*% m_prod }
                liks[anc, ] <- m_prod / comp[anc]
            }
            lk[counter] <- rowSums(log(matrix(comp[-TIPS],1)))
        }
	}

	if(is_inv){
		if(is_inv_site == T){
			prop <- c(inv, 1-inv)
            prop <- c(inv, (1-inv)/2, (1-inv)/2) # this is for G{2}+I
			log_lk <- -log(sum(prop * exp(lk)))
		} else {
			lk[is.nan(lk)] <- 0
			prop <- c(0, 1-inv)
            prop <- c(0, (1-inv)/2, (1-inv)/2)
            log_lk <- -log(sum(weight * exp(lk))) # best-practice
			#log_lk <- -log(sum(prop * exp(lk)))
		}
	}else{
        # EX2
        #freq <- c(0.448, 0.552)
        #freq <- c(0.293, 0.707)
        # EX3
        #freq <- c(0.415, 0.389, 0.196)
        #freq <- c(0, 0.985, 0.015)
        #weight <- outer(prop, freq, "*"); weight <- c(t(weight))
		#log_lk <- -log(sum(1/length(rs) * exp(lk))) # best-practice
		log_lk <- -log(sum(weight * exp(lk))) # best-practice
	}

    #print(log_lk)
	return (ifelse(is.na(log_lk), Inf, log_lk))
}


format_do_phylo_log_lk <- function(i, param){
    x <- v[[i]]
    #Qs <- list(qs[[i]])
    #Pis <- list(pis[[i]])

	liks <- matrix(0, nb.tip + nb.node, nl)
	liks[cbind(TIPS, x[[1]])] <- 1
	is_inv_site <- ifelse(all(x[[1]] == x[[1]][1]), T, F)
	environment(do_phylo_log_lk) <- environment()
	do_phylo_log_lk(param, is_inv, is_inv_site) * x[[2]]
}


sum_phylo_log_lk <- function(param){
    #cherry_liks <- calculate_lk_cherry(param)
    cherry_liks = list()
	environment(format_do_phylo_log_lk) <- environment()
    sum( unlist (mclapply(1:length(v), format_do_phylo_log_lk, param=param, mc.cores=cpu))) # PMSF
    #sum( unlist (mclapply(v, format_do_phylo_log_lk, param=param, cherry_liks=cherry_liks, mc.cores=cpu)))
}


calculate_lk_cherry <- function(param){
    cherry_liks = list()
    for (cherry in cherry_nodes){
        site_patterns <- vector()
        cherry_char <- as.character(cherry)
        cherry_liks[[as.character(cherry)]] <- list()
        for (i in 1:4) {
            for (j in 1:4) {
                vector1 <- numeric(4); vector1[i] <- 1
                vector2 <- numeric(4); vector2[j] <- 1
                lik_tmp <- matrix(c(vector1, vector2), nrow = 2)
                m<-((length(phy$tip.label)+phy$Nnode) - cherry + 1)*2;
                con_lk <- E(Q*param[m-1]) %*% vector1 * E(Q*param[m]) %*% vector2 # here 1 and 2 are inaccurate
                cherry_liks[[cherry_char]] <- c( cherry_liks[[cherry_char]], list(con_lk) )
                site_patterns <- append(site_patterns, paste(c(vector1,vector2),collapse=''))
            }
        }
        names(cherry_liks[[cherry_char]]) <- site_patterns
    }
    return(cherry_liks)
}


######################################
phy <- NULL
s <- NULL
is_inv <- F
cpu <- 1
type <- 'DNA'

is_pmsf <- F
julia_outdir <- NULL
is_force <- FALSE


######################################
spec <- matrix(c(
	'tree', 't', 2, 'character',
	'format', 'F', 2, 'character',
	'traits', 's', 2, 'character',
	'inv', 'i', 0, 'logical',
	'cpu', 'n', 2, 'integer',
    'type', '', 2, 'character',
    'pmsf', 'p', 0, 'logical',
    'julia_outdir', 'j', 2, 'character',
    'force', 'f', 0, 'logical'
	), ncol=4, byrow=T
)

opts <- getopt(spec)

if(!is.null(opts$tree)){
	phy <- read.tree(opts$tree)
    phy <- unroot(phy) # 2024.10.26
}

if(!is.null(opts$format)){
	format <- opts$format
}
if(!is.null(opts$traits)){
	s <- read.fasta(opts$traits)
	#s <- as.list(read.alignment(opts$traits, format=format)) # unsuccessful, getLength() not working on this class
}

if(!is.null(opts$inv)){
	is_inv <- T
}
if(!is.null(opts$cpu)){
	cpu <- opts$cpu
}
if(!is.null(opts$type)){
	type <- opts$type
}
if(! is.null(opts$pmsf)){
    is_pmsf <- T
}
if(!is.null(opts$julia_outdir)){
    julia_outdir <- opts$julia_outdir
}
if(!is.null(opts$force)){
    is_force <- opts$force
}


######################################
# mkdir_with_force
if(!is.null(julia_outdir)){
    if(!dir.exists(julia_outdir)){
        dir.create(julia_outdir, recursive=T)
    } else if(is_force){
        unlink(julia_outdir, recursive=T)
        dir.create(julia_outdir, recursive=T)
    } else{
        stop(paste("julia_outdir", julia_outdir, "already exists! Exiting ......", sep=" "))
    }
}

stopifnot(class(phy) != "Phylo")


######################################
#v <- read_traits(d)
if(var(getLength(s)) != 0){
	stop(paste("seqfile", "seq of diff lengths"))
}
seq <- getSequence(s) # for DNA

v <- list()
for(i in 1:getLength(s)[1]){
	v[[i]] <- sapply(s, function(x) x[i])
	v[[i]] <- v[[i]][match(phy$tip.label, names(v[[i]]))] # make them the same order
}

nl <- nlevels(factor(unlist(v, use.names=F)))


######################################
#print(v); q()
if(type == "DNA"){
    v <- lapply(v, rename_trait)
    cat("Alignment length:", length(v), "\n")
    v <- lapply(v, function(x) list(x, 1))
    site <- seq_along(v)
    cat("No. of site patterns:", length(v), "\n")
} else if(type == "AA"){
    v <- lapply(v, function(x) unname(AA_list[toupper(x)]) )
    for (i in seq_along(v)) {
        v[[i]] <- replace_not_in_AAs(v[[i]])
    }
    v <- lapply(v, function(x) list(x,1))
    dedup_rv <- deduplicate_pattern(v)
    v <- dedup_rv$v2; site <- dedup_rv$site
} else{
    stop("unknown seq type! Exiting ......")
}

#print(site); q()


######################################
nb.tip <- length(phy$tip.label)
nb.node <- phy$Nnode

TIPS <- 1:nb.tip
e1 <- phy$edge[, 1]
e2 <- phy$edge[, 2]
#print(e1)
#print(e2)
#cat(paste("nb.tip", nb.tip, "\n"))
#cat(paste("nb.node", nb.node, "\n"))

if(F){
    if(type == "DNA"){
        Q <- matrix(rep(p, 4), nl, nl)
        diag(Q) <- 0
        Q <- apply(Q, 1, function(x){x/sum(x)})
        diag(Q) <- -rowSums(Q)
        print(Q)
    } else if (type == "AA"){
        source(paste0(DIR, 'lib/read_AA_model.R'), chdir=T)
        Q_P_list <- get_LG()
        Qs <- Q_P_list$Qs
        Pis <- Q_P_list$Pis
    }
}


if(is_pmsf){
    df <- read.table('2.sitefreq')
    p <- rep(1, 20)
    Q <- matrix(rep(p, 20), nl, nl); diag(Q) <- 0; Q <- apply(Q, 1, function(x){x/sum(x)}); diag(Q) <- -rowSums(Q)
    #print(Q);
    qs <- list(); pis <- list()
    for(i in 1:nrow(df)){
        qs[[i]] <- get_new_Q(Q, unname(df[i,]))
        pis[[i]] <- as.vector(t(unname(df[i,])[-1]))
    }
    #print(pis[[1]])
    #print(sum(pis[[1]] * diag(qs[[1]])))
}

#U<-eigen(Q)$vectors; Lambda <- eigen(Q)$values; U_inv <- solve(U)

all_children <- Children(phy, (1 + nb.tip):(nb.node + nb.tip) ) #GPT
#print(all_children); q()

cherry_nodes <- find_cherry_nodes(phy)

#sum_phylo_log_lk <- compiler::cmpfun(sum_phylo_log_lk)

output_julia(julia_outdir, all_children, v)

q()


