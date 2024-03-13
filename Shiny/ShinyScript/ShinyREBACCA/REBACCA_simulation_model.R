#####################
# Generative models #
#####################

##########################
# log ratio normal model #
##########################
gm_LRN = function(Nk, N_ind=100, Size=3000, K=10^6, Sigma){
	if(length(Size)==1){
		Size = Size*rep(1,N_ind)
	}
	if(length(K)==1){
		K = rep(K, Nk)
	}
	
    # 1. convert to log-ratio basis
	K.LRN = K[1:(Nk-1)]/K[Nk]
	MVN_mu = log(K.LRN)
	# 2. multi-normal sampling
	A = diag(Nk)
	A[,Nk] = -1
	A = A[1:(Nk-1),]
	Cov.LogProp = Sigma

	MVN_sigma = A%*%Cov.LogProp%*%t(A)  #covert cov of log prop to cov of log ratio
	MVN_sample = mvrnorm(N_ind, mu=MVN_mu, Sigma=MVN_sigma)

	# 3. calculate logistic transformation
	Theta = exp(MVN_sample)
	Theta = apply(Theta, 1, function(h)h/(1+sum(h)))
	Theta = rbind(Theta, apply(Theta, 2, function(h)(1-sum(h))))

	# 4. multinomial samples
	LRM_counts = matrix(0,nrow=Nk,ncol=N_ind)
	for(i in 1:N_ind){
		LRM_counts[,i] = as.numeric(rmultinom(1, Size[i], Theta[,i]))
	}
	LRM_prop = LRM_counts%*%diag(1/Size)
	LRM_prop
}

############################ 
# Poisson log normal model #
############################
gm_LNP = function(Nk, N_ind=100, Size=3000, K=10^6, Sigma){
	if(length(Size)==1){
		Size = Size*rep(1,N_ind)
	}
	if(length(K)==1){
		K = rep(K, Nk)
	}
	
	# 1. multi-normal sampling
	MVN_mu = rep(0, Nk)
	MVN_sigma = Sigma
	#--- let component 3 and 5 be correlated
	MVN_sample = mvrnorm(N_ind, mu=MVN_mu, Sigma=MVN_sigma)
		
	# 2. generate individual propotions using poisson distribution, accounting individual variations
	pois.prop = matrix(0,nrow = Nk, ncol = N_ind)
	for (i in 1:N_ind){
		for (j in 1:Nk){
			Temp_K = K[j]*exp(MVN_sample[i,j])
			Temp_draw = rpois(1, Temp_K) 
			if(is.na(Temp_draw)){
				# redraw with normal approximation
				pois.prop[j,i] = round(rnorm(1, mean=Temp_K, sd=sqrt(Temp_K)))
			}else{
				pois.prop[j,i] = Temp_draw		
			}
		}
	}
	pois.prop = apply(pois.prop, 2, function(h)h/sum(h))

	# 3. generate samples, accounting measurement variations
	X_c = matrix(0, nrow = Nk, ncol = N_ind)
	for(i in 1:N_ind){
		X_c[,i] = rmultinom(1, Size[i], pois.prop[,i])
	}
	X_prop = X_c%*%diag(1/Size)
	X_prop
}

##############################
# Dirichlet log normal model #
##############################
gm_LND = function(Nk, N_ind=100, Size=3000, K=10^6, Sigma){
	if(length(Size)==1){
		Size = Size*rep(1,N_ind)
	}
	if(length(K)==1){
		K = rep(K, Nk)
	}	
	# 1. log-normal
	MVN_mu = rep(0, Nk)
	MVN_sigma = Sigma

	MVN_sample = mvrnorm(N_ind, mu=MVN_mu, Sigma=MVN_sigma)
	LND_sample = apply(exp(MVN_sample),1,function(h)h*K)
	
	# 2b. DM samples
	DM_samples = matrix(0, nrow = length(K), ncol = N_ind)
	for(i in 1:N_ind){
		DM_samples[,i] = DM_sampler(LND_sample[,i], Size[i])
	}
	X_prop = DM_samples%*%diag(1/Size)
	X_prop
}

# Generate Dirichlet-multinomial samples #
DM_sampler = function(Alpha, Size){
	Sample.DIR = rdirichlet(1, Alpha)
	Sample.MULTI = rmultinom(1, Size, Sample.DIR)
	Sample.MULTI
}