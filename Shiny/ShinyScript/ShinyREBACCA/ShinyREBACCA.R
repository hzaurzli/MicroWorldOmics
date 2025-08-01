library(shiny)
library(igraph)


GetMatSeq_A = function(N_seq){
  Seq_mat = matrix(0,nrow=N_seq,ncol=1/2*N_seq*(N_seq-1))
  for(i in 1:N_seq){
    Temp1=matrix(0,nrow=N_seq,ncol=N_seq)
    Temp1[i,]=1
    diag(Temp1)=0
    Temp2 = Temp1+t(Temp1)
    S_c = numeric(0)
    for(j in 1:(N_seq-1)){
      S_c = c(S_c, Temp2[1:(N_seq-j), N_seq-j+1])
    }
    Seq_mat[i, ] = S_c
  }
  Seq_mat
}

#fill into matrix sequence
# a matrix of values storing the index of a sequence
fill_seq2mat = function(N_comb){
  K = N_comb-1
  S_k = 1
  rr = matrix(0, nrow=N_comb, ncol=N_comb)
  while(K>=1){
    rr[1:K,K+1] = S_k:(K + S_k - 1)
    S_k = S_k + K
    K = K-1
  }
  rr
}

# exhange index i with j
# resulting matrix of values storing the new index of sequence
readseq_rotation = function(Temp_mat, i, j){
  Temp_rr = Temp_mat
  Temp_rr[ ,i] = Temp_mat[ ,j]
  Temp_rr[ ,j] = Temp_mat[ ,i]
  Temp = Temp_rr
  Temp_rr[i, ] = Temp[j, ]
  Temp_rr[j, ] = Temp[i, ]
  Temp = Temp_rr
  Temp = t(Temp)
  Temp[lower.tri(Temp)]=0
  Temp_rr = Temp + Temp_rr
  Temp_rr[lower.tri(Temp_rr)]=0
  Temp_rr
}

# read values from last column of matrix, top to diagonal
readseq_mat = function(SeqMat){
  N_comb = nrow(SeqMat)
  K = N_comb-1
  r_seq = numeric(0)
  while(K>=1){
    r_seq = c(r_seq, SeqMat[1:K,K+1])
    K = K-1
  }
  r_seq
}

getRotation_G = function(N_A){
  Ncol_mat = N_A*(N_A-1)/2
  Temp_rr = matrix(0,nrow=Ncol_mat,ncol=Ncol_mat)	
  #construct the first D-1 rows of equations with pairs {(D,1),...,(D,D-1)}
  Temp = matrix(rep(-1/(N_A-2),(N_A-1)*(N_A-1)),nrow=N_A-1)
  diag(Temp)=0
  Temp1 = diag(N_A-1)+Temp
  ddA = 2/((N_A-3)*(N_A-2))
  Temp2 = matrix(ddA, nrow=N_A-1, ncol=(N_A-1)*(N_A-2)/2)
  Temp3 = -(N_A-1)/((N_A-2)*(N_A-3))*GetMatSeq_A(N_A-1)
  Temp4 = cbind(Temp1, Temp2+Temp3)
  i = N_A	
  sum_i = 0
  # exchange index in a pair {(i,1),...,(i,D-1)} with {(D,1),...,(D,D-1)}
  # which is equivalent to rotate columns of matrix of Temp4
  while(i > 1){
    if(i == N_A){
      Rot_seq = 1:Ncol_mat
    }else{
      Rot_seq = readseq_mat( readseq_rotation(fill_seq2mat(N_A), i, N_A) ) 
    }
    # combine the rows and 
    # choose only the first i-1 pairs {(i,1),...,(i,i-1)}
    Temp_rr[(sum_i+1):(sum_i+i-1), ] = Temp4[1:(i-1), Rot_seq]
    sum_i = sum_i + i - 1
    i = i - 1
  }
  Temp_rr
}

# b) calculate W, the left hand side of equation (7) #
MatSumDiag_log = function(Datalog_prop){
  N_comp = nrow(Datalog_prop)
  S = matrix(0, nrow=N_comp, ncol=N_comp)
  for(i in 1:(nrow(Datalog_prop)-1)){
    for(j in i:nrow(Datalog_prop)){
      LR.temp = Datalog_prop[i,]-Datalog_prop[j,]
      S[i,j] = var(LR.temp, na.rm=T)
    }
  }
  S = S + t(S)
  diag(S) = NA
  S
}

compute_W_log = function(Data){
  N_data = nrow(Data)
  
  #compute W
  K = N_data
  S_mat = MatSumDiag_log(Data)
  S_total = sum(S_mat, na.rm=T)
  S_B = S_total
  N_A = N_data-1
  N_B = N_data
  N_C = N_data-2
  N_D = N_data-1
  W = rep(0, N_data*(N_data-1)/2)
  ID2 = 0
  
  #vector W in the order of pair combinations
  # for K components
  # {(K,1),(K,2),(K,3),...,(K,K-1),(K-1,1),...,(K-1,K-2),...,(2,1)}	
  while(K>1){
    S_A = S_B - 2*sum(S_mat[K, ], na.rm=T)
    Sigma_z1 = S_B/(2*(N_B-1))-S_A/(2*(N_A-1))
    Delta_temp = numeric(0)
    # K replaced by 10
    T_Seq = 1:N_data
    T_Seq[K] = N_data
    T_Seq = T_Seq[-N_data]
    for(i in T_Seq[1:(K-1)]){
      S_C = S_A - sum(S_mat[i,-K], na.rm=T)*2
      S_D = S_C + sum(S_mat[K,-c(K,i)], na.rm=T)*2
      Sigma_z2 = S_D/(2*(N_D-1))-S_C/(2*(N_C-1))
      Delta_temp = c(Delta_temp, 1/2*N_A*(Sigma_z2-Sigma_z1))
    }
    
    ID1 = 1 + ID2
    ID2 = ID1 + K - 2
    W[ID1:ID2] = Delta_temp
    K = K - 1
  }
  W	
}	



Delta2Cov = function(N_comb, Delta){
  K = N_comb-1
  S_k = 1
  Cov_basis = matrix(0, nrow=N_comb, ncol=N_comb)
  while(K>=1){
    Cov_basis[1:K,K+1] = Delta[S_k:(S_k+K-1)]
    S_k = S_k+K
    K = K-1
  }
  Cov_basis + t(Cov_basis)
}


rebacca = function(x, nbootstrap=50, N.cores = 1){
  dimx = dim(x)
  n = dimx[2]
  p = dimx[1]
  
  #replace zeros with minimum/10
  x_min = min(x[x!=0])
  x[x==0] = x_min/10
  
  # if data is not proportions, normalize the data
  if(max(x)>1){
    x = x%*%diag(1/apply(x, 2, sum))
  }
  
  #log transform
  x = log(x)
  
  cat("calculating rotation matrix\n")	
  tempG = getRotation_G(p)
  
  #pre-calculate LASSO path using full data
  tempW = compute_W_log(x)
  pre.lasso = tryCatch(glmnet(tempG, tempW, alpha=1, standardize = F, intercept=F), error=function(
    e){cat("obtaining pre-LASSO fails!\n")})
  Lambda = pre.lasso$lambda		
  
  #choose lambda path to reveal proper amount of nonzero solutions
  # maximum df = 1/2 of total variables
  q_max = 1/2*p*(p-1)/2
  
  Temp = which(pre.lasso$df<q_max)
  Lambda_max_id = Temp[length(Temp)]
  Lambda = Lambda[1:Lambda_max_id]
  
  # initial values
  freq = matrix(0, p*(p-1)/2, length(Lambda))
  qval = 0
  S_Mat = numeric(0)
  
  # parallel calculation
  cat("calculating LASSO with parallel bootstrapping\n")
  mc = getOption("mc.cores", N.cores)
  res = mclapply(seq_len(N.cores), select_parallel, Data=x, RotationG=tempG, Lambda=Lambda, N_boots=nbootstrap, N_cores=mc, mc.cores=mc)
  
  # merge all paralleled results 
  # add frequencies up 
  for (i in 1:length(res)){
    freq = freq + res[[i]][[1]]
    qval = qval + res[[i]][[2]]
  }
  
  # normalize frequence in [0,1]
  freq = freq/(2*nbootstrap)
  qval = qval/(2*nbootstrap)
  
  # choose the maximum selection probability for each variable
  S_score = apply(freq, 1, max)
  
  # convert to matrix with significance
  RR = Delta2Cov(p, S_score)
  diag(RR) = NA
  
  list(Stability = RR, q = qval)
}

# multi-core running code
select_parallel = function(cc, Data, RotationG, Lambda, N_boots, N_cores){
  bin_size = floor(N_boots/N_cores)
  seq_start = (cc-1)*bin_size + 1
  seq_end = ifelse(cc==N_cores, N_boots, cc*bin_size)
  n = ncol(Data)
  p = nrow(Data)
  # global variable tempG and Lambda
  freq = matrix(0, p*(p-1)/2, length(Lambda))
  q_temp = 0
  j = 1
  
  for (i in seq_start:seq_end){
    perm = sample(n, replace=F)
    
    #use split-sample resampling method, stabiity selection		
    #calculate W vector on each set of subsamples
    tempW1 = compute_W_log(Data[, perm[1:floor(n/2)]])
    tempW2 = compute_W_log(Data[, perm[floor(n/2):n]])
    
    r1 = tryCatch(glmnet(RotationG, tempW1, alpha=1, standardize = F, intercept=F, lambda = Lambda), error=function(
    e){cat("warning: fails to run LASSO")})
    r2 = tryCatch(glmnet(RotationG, tempW2, alpha=1, standardize = F, intercept=F, lambda = Lambda), error=function(
    e){cat("warning: fails to run LASSO")})
    
    # extract beta from lars (need transpose) or glm
    beta1 = abs(sign(r1$beta))
    beta2 = abs(sign(r2$beta))
    Temp1 = apply(beta1, 1, function(h)ifelse(sum(h)==0, 0, 1))
    Temp2 = apply(beta2, 1, function(h)ifelse(sum(h)==0, 0, 1))	
    q_temp = q_temp + length(which(Temp1 != 0)) + length(which(Temp2 != 0))
    freq = freq + beta1 + beta2
    #SelectMat[ ,j] = Temp1
    j = j + 1
  }
  list(freq, q_temp)
}


r.TailProbs <- function(eta, B, r) {
  # TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
  # We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
  MAXa <- 100000
  MINa <- 0.00001
  s <- -1/r
  etaB <- eta * B
  k_start <- (ceiling(2 * etaB) + 1)
  output <- rep(1, B)
  if(k_start > B) return(output)
  
  a_vec <- rep(MAXa,B)
  
  Find.a <- function(prev_a) uniroot(Calc.a, lower = MINa, upper = prev_a, tol = .Machine$double.eps^0.75)$root
  Calc.a <- function(a) {
    denom <- sum((a + 0:k)^(-s))
    num <- sum((0:k) * (a + 0:k)^(-s))
    num / denom - etaB
  }
  
  for(k in k_start:B) a_vec[k] <- Find.a(a_vec[k-1])
  
  OptimInt <- function(a) {
    num <- (k + 1 - etaB) * sum((a + 0:(t-1))^(-s))
    denom <- sum((k + 1 - (0:k)) * (a + 0:k)^(-s))
    1 - num / denom
  }
  # NB this function makes use of several gloabl variables
  
  prev_k <- k_start
  for(t in k_start:B) {
    cur_optim <- rep(0, B)
    cur_optim[B] <- OptimInt(a_vec[B])
    if (prev_k <= (B-1)) {
      for (k in prev_k:(B-1))
        cur_optim[k] <- optimize(f=OptimInt,lower = a_vec[k+1], upper = a_vec[k], maximum  = TRUE)$objective
    }
    output[t] <- max(cur_optim)
    prev_k <- which.max(cur_optim)
  }
  return(output)
}

minD <- function(theta, B, r = c(-1/2, -1/4)) {
  pmin(c(rep(1, B), r.TailProbs(theta^2, B, r[1])), r.TailProbs(theta, 2*B, r[2]))
}


stability_cutoff = function(Stab_score, q_val, B=50, FWER=0.05){
  p = nrow(Stab_score)
  N_p = p*(p-1)/2
  Pval_ref = minD(q_val/N_p, B)
  
  #largest one that smaller than FWER=0.01
  Temp_select = which(Pval_ref < FWER)[1]
  Temp_select/(B*2)
}

#convert stability score to adjacency matirx
sscore2adjmatrix = function(S_mat, alpha){
  S_mat[S_mat < alpha] = NA
  S_mat[!is.na(S_mat)] = 1
  S_mat[is.na(S_mat)] = 0
  as.matrix(S_mat)
}

#convert cov matrix to variables
Cov2Delta = function(Cov_basis){
  N_sp = nrow(Cov_basis)
  N_p = N_sp*(N_sp-1)/2
  Delta = rep(0, N_p)
  K = N_sp - 1
  S_k = 1
  while(K>=1){
    Delta[S_k:(S_k+K-1)] = Cov_basis[1:K,K+1]
    S_k = S_k+K
    K = K-1
  }
  Delta
}


rebacca_adjm2corr = function(Data, Adj_Mat){
  x = as.matrix(Data)
  dimx <- dim(x)
  n <- dimx[2]
  p <- dimx[1]
  
  #replace zeros with minimum/10
  x_min = min(x[x!=0])
  x[x==0] = x_min/10
  
  # if data is not proportions, normalize the data
  if(max(x)>1){
    x = x%*%diag(1/apply(x, 2, sum))
  }
  
  #log transform
  x = log(x)
  
  cat("calculating rotation matrix\n")	
  TempG = getRotation_G(p)
  TempW = compute_W_log(x)
  
  cat("calculating ls fit\n")
  Temp_b0 = Cov2Delta(Adj_Mat)
  uncons_ind = which(Temp_b0==1)
  
  #remove columns in G that correspond the variables with zeros
  TempG2 = TempG[, uncons_ind]
  
  #use least square fit for nonzero variables
  Temp = lm.fit(y=TempW, x=TempG2)
  
  # set nonzero variables equals lm fit estimates
  Estimate_lmfit = Temp_b0
  Estimate_lmfit[uncons_ind] = Temp$coef
  
  #----calculate correlations----#
  cat("calculating correlation\n")
  Result_cov = Delta2Cov(p, Estimate_lmfit)
  
  #calculate variance
  S_mat = MatSumDiag_log(x)
  S_total = sum(S_mat, na.rm=T)
  
  TotalSum = (S_total + 2*sum(Result_cov))/(2*(p-1))
  Var_basis = rep(0, p)
  
  for(i in 1:p){
    Data_ratio =t(apply(x, 1, function(h)(h-x[i,])))
    LR_temp = diag(cov(t(Data_ratio), use="complete.obs"))
    Var_basis[i] = (sum(LR_temp) - TotalSum + 2*sum(Result_cov[i,]))/(p-2)
  }
  diag(Result_cov) = Var_basis
  list(cov = Result_cov, corr = cov2cor(Result_cov))
}


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


ui <- tagList(
  fluidPage(
    titlePanel("Network analysis based REBACCA"),
    sidebarLayout(
      sidebarPanel(
        # uiOutput 做上传文件的 ui, 对应后面的 output$file1
        uiOutput('file1'),
        actionButton('reset', 'RESET'),
        actionButton('start', 'START'),
        hr(),
        downloadButton("downloadData", "Download Table"),
        hr(),
        h5('Developer:'),
        h6('Small runze (shiny app)'),
        br(),
        h5('Github: '),
        h6('https://github.com/hzaurzli (Small runze)'),
      ),
      mainPanel(
        h4("Relationship matrix"),
        shinycssloaders::withSpinner(
          dataTableOutput("table")
        ),
        hr(),
        plotOutput(outputId = "detectfig")
        
      )
    )
  )
)



server <- function(input, output, session) {
  options(shiny.maxRequestSize=1024*1024*1024^2)
  
  values <- reactiveValues(
    file = NULL
  )
  
  
  otu_dat <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,header = T,row.names = 1)
  })
  
  
  
  # observeEvent(input$reset), 代表点击 RESET 时触发的动作,此时重新渲染 fileInput 的 ui
  observeEvent(input$reset, {
    values$file <- NULL
    output$file1 <- renderUI({
      fileInput("file1", "Step 1: Choose otu abundance matrices",
                accept=c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
      )
    })
  }, ignoreNULL = F)
  
  
  output$table <- renderDataTable({
    data = data.frame(Item = c('No data'))
  }, options = list(pageLength = 1, searching = FALSE, paging = FALSE))
  
  
  observeEvent(input$start, {
    
    otu_dat = otu_dat()
    
    if(is.null(otu_dat)){
      warning("Please upload files!")
    }
    else{
      output$table <- renderDataTable({
        require(MCMCpack)
        require(glmnet)
        require(parallel)
        
        x = as.matrix(otu_dat)
        x.rslt = rebacca(x, nbootstrap=100, N.cores=1)
        
        # estimate correlation
        tau = stability_cutoff(x.rslt$Stability, x.rslt$q, B=50, FWER=0.05)
        x.adj = sscore2adjmatrix(x.rslt$Stability, tau)
        x.est = rebacca_adjm2corr(x, x.adj)
        
        myAdjacencyMatrix = x.est$corr
        
        g  <- graph.adjacency(myAdjacencyMatrix,weighted=TRUE)
        dat <<- get.data.frame(g)
      },options = list(pageLength = 10))
    }
  })
  
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(dat,file,row.names = T,quote = F)
    }
  )
}


# Run the application 
app = shinyApp(ui = ui, server = server)
runApp(app, port = 50737)
