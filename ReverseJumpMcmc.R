library(nimble)
library(coda)
library(bayestestR)
############################################################
###########             Load data               ############
############################################################

TD<-read.csv(file="treatment_death_data.csv",header=T)

#DeathData<-read.csv(file="../../DeathData_07-17.csv",header=T)
#y.D = apply(DeathData[,3:9],1,sum)
y.D = TD$Death
y.T = TD$Treatment_u21+TD$Treatment_o21
Pop = TD$TotalPopulation

### identify what counties are censored
cenIndc = 1*(is.na(TD$Treatment_u21) & !is.na(TD$Treatment_o21))
cenInda = 1*(!is.na(TD$Treatment_u21) & is.na(TD$Treatment_o21))
cenIndb = 1*(is.na(TD$Treatment_u21) & is.na(TD$Treatment_o21))
cenlb = rep(NA,length(y.T))
cenub = rep(NA,length(y.T))
cenlb[which(cenIndc==1)]=TD$Treatment_o21[which(cenIndc==1)]+1
cenub[which(cenIndc==1)]=TD$Treatment_o21[which(cenIndc==1)]+9
cenlb[which(cenInda==1)]=TD$Treatment_u21[which(cenInda==1)]+1
cenub[which(cenInda==1)]=TD$Treatment_u21[which(cenInda==1)]+9
cenlb[which(cenIndb==1)]=2
cenub[which(cenIndb==1)]=18
cenlb[which(is.na(cenlb))]=TD$TotalPopulation[which(is.na(cenlb))]
cenub[which(is.na(cenub))]=TD$TotalPopulation[which(is.na(cenub))]



### simply imputing values for the counties that were missing
### can add more formal methods later
y.Timp = y.T
y.Timp[which(cenIndc==1)] = TD$Treatment_o21[which(cenIndc==1)]+5
y.Timp[which(cenInda==1)] = TD$Treatment_u21[which(cenInda==1)]+5
y.Timp[which(cenIndb==1)] = 10


W<-read.csv(file = "OhioAdjacency.csv",header = FALSE)
n<-nrow(W)

ET = rep(0,n)
ED = rep(0,n)
for(i in 1:n){
  ET[i] =  sum(y.Timp)/sum(Pop)*Pop[i]
  ED[i] =  sum(y.D)/sum(Pop)*Pop[i]
}

### construct design matrix of covariates
X.T = cbind(log(100-TD$PercentWhite),log(TD$MedianHouseholdIncome), TD$Disability18_64,TD$PublicInsurance,TD$HPSA)
for(i in 1:length(X.T[1,])){
  X.T[,i] = (X.T[,i]-mean(X.T[,i]))/sd(X.T[,i])
}
X.D = cbind(log(TD$MedianHouseholdIncome),TD$MedianAge,TD$FemaleHouse)
for(i in 1:length(X.D[1,])){
  X.D[,i] = (X.D[,i]-mean(X.D[,i]))/sd(X.D[,i])
}
X.f = cbind(log(TD$Bachelors),TD$HIDTA,TD$BelowPoverty,1*(TD$Interstates>0))
for(i in 1:length(X.f[1,])){
  X.f[,i] = (X.f[,i]-mean(X.f[,i]))/sd(X.f[,i])
}
bp = length(X.f[1,])
bp.T = length(X.T[1,])
bp.D = length(X.D[1,])


############################################################
###########           Set up for NIMBLE           ############
############################################################


adj<-NULL
for(j in 1:n){
  adj<-c(adj,which(W[j,]==1))
}
num = colSums(W)
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*adj

cen=cenInda+cenIndc+cenIndb
bd=cbind(cenlb,cenub)

model_code=nimbleCode({
  
  ## create a z container to mark importance of variables
  for(i in 1:bp) {
    z[i] ~ dbern(psi[i]) ## indicator variable for each coefficient
    beta[i] ~ dflat()
    #psi[i] ~ dunif(0,1)
    psi[i] <- 0.5
    zbeta[i] <- z[i] * beta[i]  ## indicator * beta
  }
  
  for(i in 1:bp.T) {
    z.T[i] ~ dbern(psi.T[i]) ## indicator variable for each coefficient
    beta.T[i] ~ dflat()
    #psi.T[i] ~ dunif(0,1)
    psi.T[i] <- 0.5
    zbeta.T[i] <- z.T[i] * beta.T[i]  ## indicator * beta
  }
  
  for(i in 1:bp.D) {
    z.D[i] ~ dbern(psi.D[i]) ## indicator variable for each coefficient
    beta.D[i] ~ dflat()
    #psi.D[i] ~ dunif(0,1)
    psi.D[i] <- 0.5
    zbeta.D[i] <- z.D[i] * beta.D[i]  ## indicator * beta
  }
  
  for (i in 1:n){
    y.D[i] ~ dpois(ED[i]*lambda1[i])
    y.T[i] ~ dpois(ET[i]*lambda2[i])
    cen[i] ~ dinterval(y.T[i],bd[i,1:2])
    
    ## regression equations
    log(lambda1[i])<- beta0.D+inprod(X.D[i,1:bp.D],zbeta.D[1:bp.D])+(aa+alpha[i])*(U[i]+mu.U[i])+ V1[i]
    log(lambda2[i])<- beta0.T+inprod(X.T[i,1:bp.T],zbeta.T[1:bp.T])+(U[i]+mu.U[i])+  V2[i]
    mu.U[i] <- inprod(X.f[i,1:bp],zbeta[1:bp])
    
    V1[i] ~ dnorm(0,tau.V1) 
    V2[i] ~ dnorm(0,tau.V2)
  }
  
  # ICAR prior for the spatial random effects and loadings alpha
  U[1:n] ~ dcar_normal(adj[], weights[], num[], tau.U,zero_mean=1)
  alpha[1:n] ~ dcar_normal(adj[], weights[], num[], tau.alpha,zero_mean=1)
  
  #Priors
  beta0.T ~ dflat()
  beta0.D ~ dflat()

  tau.V1 ~ dgamma(0.5,0.5)
  tau.V2 ~ dgamma(0.5,0.5)
  aa ~ dnorm(1,.5)
  tau.U ~ dgamma(0.5,0.5)
  tau.alpha ~ dgamma(0.5,0.5)
})


############################################################
#########              Call NIMBLE                 ###########
############################################################

TD_constants=list(n=n,bp = length(X.f[1,]),bp.T = length(X.T[1,]),bp.D=length(X.D[1,]),num=num,adj=adj,weights=weights,bd=bd)

TD_data=list(y.D=y.D,y.T=y.T,ET=ET,ED=ED,X.f=X.f,X.D=X.D,X.T=X.T,cen=cen)


#glm0 = glm(TDbyyear$AllDeaths~1,offset=log(TDbyyear$TotalPopulation))

TD_inits=list(beta0.T=0,beta.T=rep(0,length(X.T[1,])),beta0.D=0,beta.D=rep(0,length(X.D[1,])),beta = rep(0,length(X.f[1,])),U=rep(0,n),
              V1=rep(0,n),V2=rep(0,n),alpha=rep(0,n),aa=1,
              tau.V1=1,tau.V2=1,tau.U=1,tau.alpha=1,y.T=y.Timp,
              z = rep(1,bp), psi = rep(.9,bp),
              z.T=rep(1, bp.T), psi.T = rep(.9,bp.T), 
              z.D=rep(1, bp.D), psi.D = rep(.9,bp.D))


# Build the model.
TD_model <- nimbleModel(model_code, TD_constants,TD_data,TD_inits)
TD_compiled_model <- compileNimble(TD_model,resetFunctions = TRUE)

# Set up samplers.
TD_mcmc_conf <- configureMCMC(TD_model,monitors=c("beta0.T","beta0.D","zbeta.T","zbeta.D","zbeta","tau.V1","tau.V2","tau.U","alpha","tau.alpha",
                                                  "U","aa","V1","V2","z", "z.T", "z.D","lambda2","lambda1","ET","ED","y.T"),useConjugacy = TRUE)

#TD_mcmc_conf <- configureMCMC(TD_model,useConjugacy = TRUE)


TD_mcmc<-buildMCMC(TD_mcmc_conf, enableWAIC = TRUE)
TD_compiled_mcmc<-compileNimble(TD_mcmc, project = TD_model,resetFunctions = TRUE)

# Run the model 
MM = 10000
st<-Sys.time()
TD_samples=runMCMC(TD_compiled_mcmc,inits=TD_inits,
                   nchains = 1, nburnin=MM/2,niter = MM,samplesAsCodaMCMC = TRUE,thin=5,
                   summary = FALSE, WAIC = TRUE,progressBar=TRUE) 
Sys.time()-st


TD_samples$WAIC



#save(TD_samples,file="FactorModelOutput.Rda")

samples <- as.data.frame(TD_samples)
mean(samples)

# calculate WAIC
TD_samples$calculateWAIC(nburnin = MM/2)

ci_hdi <- ci(TD_samples, method = "HDI")
# gives us the high density interval for X.f, X.D, X.T
ci_hdi[354:365, ]
ci_eti <- ci(TD_samples, method = "ETI")
# gives us the equal tailed interval for X.f, X.D, X.T
ci_eti[354:365, ]
