###########
# Example #
###########
require(MCMCpack)
require(glmnet)
require(parallel)
source("REBACCA_main.R")
source("REBACCA_simulation_model.R")

##################################################
# Part I: simulating a simple correlated network #
################################################## 

# We generate a sample of 10 components, 100 microbiome samples, sequencing depth of 3,000 counts,
# mean basis abundance 10^6 for each components
# correlation network as in simulation case 1
Nk = 10

# construct a covariance matrix according to the correlation network
c1 = c(0.3, -0.4, 0.6)
c2 = c(0.3, 0.6, -0.3)
x.Sigma = matrix(0,Nk,Nk)
for(i in 1:3){
	j= 3*i+1
	if(i==1){
		x.Sigma[i,(j-2):j] = c1
	}else{
		x.Sigma[i,(j-2):j] = c2
	}
}
x.Sigma = x.Sigma + t(x.Sigma)
diag(x.Sigma) = 1

# theoretical correlation
x.cor = cov2cor(x.Sigma)

# generate proportion data
x = gm_LRN(Nk=Nk, N_ind=100, Size=3000, K=10^6, Sigma=x.Sigma)

##############################################
# Part II: run REBACCA on the simulated data #
##############################################
x.rslt = rebacca(x, nbootstrap=50, N.cores=1)

# estimate correlation
tau = stability_cutoff(x.rslt$Stability, x.rslt$q, B=50, FWER=0.05)
x.adj = sscore2adjmatrix(x.rslt$Stability, tau)
x.est = rebacca_adjm2corr(x, x.adj)

# compare the estimated correlation with the theoretical one
x.est$corr
x.cor
