#! /usr/bin/env Rscript


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
    library(parallel)
    library(matrixStats)

    library(seqinr)
    library(phytools)
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

	if (!is.factor(x)) {
		x <- factor(x)
		nl <- nlevels(x)
		#lvls <- levels(x)
		x <- as.integer(x)
	}
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
        q()
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
	names(v[[i]]) <- names(s)
}

nl <- nlevels(factor(unlist(v, use.names=F)))


######################################
#print(v); q()
if(type == "DNA"){
    v <- lapply(v, rename_trait)
    cat("Alignment length:", length(v), "\n")
    v <- get_v_freq(v)
    cat("No. of pattern:", length(v), "\n")
} else if(type == "AA"){
    v <- lapply(v, function(x) unname(AA_list[toupper(x)]) )
    v <- lapply(v, function(x) list(x,1))
    #v <- lapply(seq, function(x){x = toupper(x); x})
} else{
    stop("unknown seq type! Exiting ......")
}

#print(v); q()


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



df <- length(phy$edge.length)
bl <- rep(0.1, df)
#bl <- c(0.1351434449,0.104783017,0.0499955681,0.0643653305,0.010691266,0.0264261792,0.1618177939,0.1461096787,0.1017798321,0.2895916068,0.2528479589,0.3269044605,0.0822965694,0.0952360872,0.0034494216,2.5259e-06,0.1153463186,0.0662615733,0.3941393636,0.0996489034,0.036545514,0.0424834086,0.0168373696,0.0243012902,0.0402086644,0.0414002208,0.0943864195,0.1042118786,0.0564096446,0.0422706142,0.0523381787,2.5259e-06,0.0596987589,0.0173806719,0.0032752347,0.0035325072,0.0562832985,2.5259e-06,0.0852478336,0.0221226782,0.1370157539,0.0966053109,0.0417246049,0.1608164433,2.5259e-06,0.0441694058,0.2333463531,0.002285576,0.0914264439,0.0598247098,0.0793496784,0.2949696846,0.0147588953,0.0063825098,0.011180306,0.0904691973,0.1311192811)
#bl <- c(1.3716642546,1.8082351105,0.887,0.339,1.556)
#bl <- c(0.2525075204,0.3598821266,0.0875467953, 0.3066884149, 0.4753269704, 0.0875467953);
#bl <- c(0.1713275668, 0.0000025259, 0.0722891959, 0.1135710432, 0.1086154433);


#bl <- c(0.0322151843, 0.0317750062, 0.0060064457, 0.0152227065, 0.0131411075, 0.0009737001, 0.0427571326, 0.0124979104, 0.0989260553, 0.1106187653, 0.1560053038, 0.0267877333, 0.0278560113, 0.0318322365, 0.0000010000, 0.3371242301, 0.2063122571)

#print ( -1 * numDeriv::hessian(sum_phylo_log_lk, bl) ); q()
#hessian <- -1 * pracma::hessian(sum_phylo_log_lk, bl); print(hessian); q()
#bl = c(0.11922693817462916, 0.07375750321775916, 0.006234746567309918, 0.0502961216239749, 0.020367838948476095, 3.859015472849246e-5, 0.20222546287553791, 0.12107420476031365, 0.31287961037130235, 0.4960120343692124, 0.3615061461320902, 0.14059651243781832, 0.10832705026318967, 0.1357459789657297, 5.543854475366512e-5, 5.321805835044186e-5, 0.1366789303462951, 0.08729706647760276, 0.6463493080864843, 0.07323815032869688, 0.015690813909650436, 0.05860574173145091, 4.7407242962270605e-5, 1.9084780864747587e-5, 0.0736377602734294, 5.4416488726763536e-5, 0.06909947109505349, 0.155412749743196, 0.030979407865874022, 0.020483959637346217, 0.030680793241475345, 0.02048295645619988, 0.08476302612482993, 4.8285572299517085e-5, 6.201060894107256e-5, 6.203637729270824e-5, 0.07380888749458947, 3.690627357307275e-5, 0.038635815673950566, 0.024624857019554695, 0.1634414121047208, 0.10951280548231543, 0.03172105996902549, 0.15812317101127146, 2.101448696132541e-5, 0.04624792197753962, 0.22090249989625255, 0.020870289945368527, 0.03184078478278123, 0.0008876157364406605, 5.3606957511628194e-5, 0.010136462566728713, 0.21893777233250766, 0.17204575670341854, 0.2591149807728129, 0.0003953739401963011, 0.6741521722715208)

#v <- v[1]
for(i in 1:1){
    #print(sum_phylo_log_lk(bl))
    print(i)
    sum_phylo_log_lk(bl)
}
q()

#q()

param <- bl


#if(is_inv) param <- append(param, 0.1)

#print(v)


#-pracma::grad(sum_phylo_log_lk, bl); q()
#opt <- nlminb(param, function(p) sum_phylo_log_lk(p), lower = 1e-6, upper = Inf, 
#	control=list(iter.max=100, eval.max=100, trace=1))

#opt = optim(par=param, sum_phylo_log_lk, method="L-BFGS-B", lower=1e-6, upper=Inf, hessian=F, control = list(maxit=100, trace=1, REPORT=1))
#print(opt$par); q()


print ( -1 * pracma::hessian(sum_phylo_log_lk, param) )

q()

-numDeriv::grad(sum_phylo_log_lk, opt$par)
-pracma::grad(sum_phylo_log_lk, c(0.092508, 0.196664, 0.019797, 0.373649, 0.569339) )
-pracma::grad(sum_phylo_log_lk, c(0.0925074480, 0.1966638668, 0.0200486974, 0.3731064023, 0.5688556173) )


