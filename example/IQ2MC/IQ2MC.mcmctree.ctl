seed = -1
seqfile = IQ2MC/IQ2MC.dummy.phy  * A dummy alignment only allow to run MCMCtree
treefile = IQ2MC/IQ2MC.rooted.nwk  * Rooted newick tree with annotated fossil/tip dates
mcmcfile = IQ2MC/IQ2MC.mcmctree.log  * MCMC log of parameters that can be examined in Tracer
outfile = IQ2MC/IQ2MC.mcmctree.out  * Output of the summerized results of MCMCtree
ckpfile = IQ2MC/IQ2MC.mcmctree.ckp  * Checkpoint file of MCMCtree
hessianfile = IQ2MC/IQ2MC.mcmctree.hessian  * File with gradient and hessian matrix


checkpoint = 1 * 0: nothing; 1 : save; 2: resume
ndata = 1 * number of partitions
seqtype = 0    * 0 : nucleotides; 1: codons; 2: AAs (not required if the approximate likelihood method is used)
usedata = 2    * 0: sampling from priors with no data; 1: exact slow likelihood; 2: approximate likelihood
clock = 3     * 1: global clock with equal rates; 2: independent rates; 3: correlated rates
RootAge = <1.0  * safe constraint on root age, used if no fossil for root.

BDparas = 1,1,0.5    * birth rate, death rate, sampling priors for sampling times
finetune = 1: 0.1  0.1  0.1  0.01 .5  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErrprint = 1  * 1: normal output; 2: verbose output

*** These parameters are used for multi-loci partitioned data (ndata > 1), see dos Reis et al .(2013)

rgene_gamma = 2 2     * alpha and beta parameter of Dirichlet-gamma prior for mean rate across loci for clock=2 or 3
sigma2_gamma = 1 10   * alpha and beta parameter of Dirichlet-gamma prior for rate variance across loci for clock=2 or 3

*** These parameters control the MCMC run

burnin = 20000
sampfreq =  100
nsample =  20000

***  Note: Total number of MCMC iterations will be burnin + (sampfreq * nsample)

***  Following parameters only needed to run MCMCtree with exact likelihood (usedata = 1), no need to change anything for approximate likelihood (usedata = 2)

model = 0      * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
alpha = 0      * 0: No rate heterogeneity across sites; otherwise: initial alpha parameter of the Gamma distribution
ncatG = 0      * Number of rate categories for the discrete Gamma distribution

cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?

kappa_gamma = 6 2      * gamma prior for kappa of the HKY model
alpha_gamma = 1 1       * alpha and beta parameter of Gamma distribution for heterogeneous rates across sites

