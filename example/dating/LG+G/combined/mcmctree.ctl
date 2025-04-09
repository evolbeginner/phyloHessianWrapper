          seed = -1
*       seqfile = tmp/seqs.phy
*      treefile = tmp/calibration_tree.txt
       seqfile = combined.phy
      treefile = species.trees
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1
       seqtype = 2
       usedata = 2 in.BV 3
         clock = 3
*      RootAge = 'B(3.20,4.51,0.025,0.01)'  * safe constraint on root age, used if no fossil for root.

         model = 0 * 0 : JC69, 1: K80, 2: F81, 3: F84, 4: HKY85
         alpha = 0.5
         ncatG = 4
   duplication = 1
     cleandata = 0    * remove sites With ambiguity data (1:yes, 0:no)?

       BDparas = 0.4 0.2 0.1
   kappa_gamma = 6 2       * gamma prior for kappa
   alpha_gamma = 1 1       * gamma prior for alpha

   rgene_gamma = 1 50 1
  sigma2_gamma = 1 10 1

      finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1): times, musigma2, rates, mixing, paras, FossilErr

         print = 2
        burnin = 1000
      sampfreq = 10
       nsample = 500
