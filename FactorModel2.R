library(nimble)
library(coda)
library(igraph)
############################################################
###########             Load data               ############
############################################################

TD<-read.csv(file="~/Desktop/summer research/Ohio data/treatment_death_data.csv",header=T)

#DeathData<-read.csv(file="../../DeathData_07-17.csv",header=T)
#y.D = apply(DeathData[,3:9],1,sum)

logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

#p<-runif(88, min = 0.25, max = 1)
#p <- logit2prob(0.5*TD$Disability18_64)

# count of death
y.D = TD$Death
# count of treatment admissions
y.T = TD$Treatment_u21+TD$Treatment_o21
# population in county i
Pop = TD$TotalPopulation

### identify what counties are censored
# standards: 
# find treatment counts < 10 for adolescents
cenIndc = 1*(is.na(TD$Treatment_u21) & !is.na(TD$Treatment_o21))
# find treatment counts < 10 for adults 
cenInda = 1*(!is.na(TD$Treatment_u21) & is.na(TD$Treatment_o21))
# find treatment counts < 10 for both adults and adolescents
cenIndb = 1*(is.na(TD$Treatment_u21) & is.na(TD$Treatment_o21))
# create an interval list to store those information together
# lower bound
cenlb = rep(NA,length(y.T)) 
# upper bound
cenub = rep(NA,length(y.T))

# we find 4 counties (only adolescent < 10)
# we find 1 counties (both < 10)
# true count lies between adolescents count and adolescents count + 9
cenlb[which(cenIndc==1)] =TD$Treatment_o21[which(cenIndc==1)]+1
cenub[which(cenIndc==1)] =TD$Treatment_o21[which(cenIndc==1)]+9
# true count lies between adult count and adult count + 9
cenlb[which(cenInda==1)] =TD$Treatment_u21[which(cenInda==1)]+1
cenub[which(cenInda==1)] =TD$Treatment_u21[which(cenInda==1)]+9
# true count lies between (adult+ adolescents) +2 and (adult+ adolescents) + 18
cenlb[which(cenIndb==1)] =2
cenub[which(cenIndb==1)] =18
cenlb[which(is.na(cenlb))] =TD$TotalPopulation[which(is.na(cenlb))]
cenub[which(is.na(cenub))] =TD$TotalPopulation[which(is.na(cenub))]



### simply imputing values for the counties that were missing
### can add more formal methods later
y.Timp = y.T
# those are initial guess for counties that have missing data
y.Timp[which(cenIndc==1)] = TD$Treatment_o21[which(cenIndc==1)]+5
y.Timp[which(cenInda==1)] = TD$Treatment_u21[which(cenInda==1)]+5
y.Timp[which(cenIndb==1)] = 10


W<-read.csv(file = "~/Desktop/summer research/Ohio data/OhioAdjacency.csv",header = FALSE)
n<-nrow(W)

# expected treatment admissions for county i
ET = rep(0,n)
# expected death for county i
ED = rep(0,n)
# compute the baseline expected values
for(i in 1:n){
  ET[i] =  sum(y.Timp)/sum(Pop)*Pop[i]
  ED[i] =  sum(y.D)/sum(Pop)*Pop[i]
}

### construct design matrix of covariates
X = cbind(log(TD$MedianHouseholdIncome),TD$Disability18_64)
X[,1] = (X[,1]-mean(X[,1]))/sd(X[,1])
X[,2] = (X[,2]-mean(X[,2]))/sd(X[,2])

############################################################
###########           Set up for NIMBLE           ############
############################################################


adj<-NULL
for(j in 1:n){
  #remeber W is the proximity
  #this helps us to find the # of neighboring counties for i
  #store the information into an adjacency matrix
  adj<-c(adj,which(W[j,]==1))
}
num = colSums(W)
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*adj

# incoporate the censored counties
cen=cenInda+cenIndc+cenIndb
bd=cbind(cenlb,cenub)

model_code=nimbleCode({
  for (i in 1:n){
    #the outcomes follow bivariate Poisson distribution
    #lambda: the relative risk of D and T in county compared to the state average
    y.D[i] ~ dpois(ED[i]*lambda1[i])
    y.T[i] ~ dpois(ET[i]*lambda2[i])
    #also include the case of the censored counties
    cen[i] ~ dinterval(y.T[i],bd[i,1:2])
    
    #joint generalized linear mixed effects model
    #beta0 are the intercepts
    #v is the common spatial factor for county i (random effects)
    #alpha are the factor loadings
    log(lambda1[i])<- beta0.D+(aa+alpha[i])*(U[i]+mu.U[i])+ V1[i]
    log(lambda2[i])<- beta0.T+(aa2+alpha2[i])*(U[i]+mu.U[i])+  V2[i]
    # compute the inner product of the matrix
    mu.U[i] <- inprod(X[i,],beta[])
    V1[i] ~ dnorm(0,tau.V1) 
    V2[i] ~ dnorm(0,tau.V2)
  }
  
  # ICAR prior for the spatial random effects v and loadings alpha
  U[1:n] ~ dcar_normal(adj[], weights[], num[], 1,zero_mean=1)
  alpha[1:n] ~ dcar_normal(adj[], weights[], num[], tau.alpha,zero_mean=1)
  alpha2[1:n] ~ dcar_normal(adj[], weights[], num[], tau.alpha2,zero_mean=1)
  
  #Priors
  #uninformative flat priors for the intercepts
  beta0.T ~ dflat()
  beta0.D ~ dflat()
  for(j in 1:bp){
    #uninformative, independent, normal, flat prior for Beta
    beta[j] ~ dflat()
  }
  #conjugate priors
  #inverse gamma for varianes
  tau.V1 ~ dgamma(0.5,0.5)
  tau.V2 ~ dgamma(0.5,0.5)
  aa ~ dnorm(1,.5)
  aa2 ~ dnorm(1,.5)
  tau.alpha ~ dgamma(0.5,0.5)
  tau.alpha2 ~ dgamma(0.5,0.5)
})


############################################################
#########              Call NIMBLE                 ###########
############################################################

TD_constants=list(n=n,bp = length(X[1,]),num=num,adj=adj,weights=weights,bd=bd)

#outcomes, expected values, explanatory matrix, 
TD_data=list(y.D=y.D,y.T=y.T,ET=ET,ED=ED,X=X,cen=cen)


#glm0 = glm(TDbyyear$AllDeaths~1,offset=log(TDbyyear$TotalPopulation))

TD_inits=list(beta0.T=0,beta0.D=0,beta = c(0,0),U=rep(0,n),
              V1=rep(0,n),V2=rep(0,n),alpha=rep(0,n),alpha2=rep(0,n),aa=1,aa2=1,
              tau.V1=1,tau.V2=1,tau.alpha=1,tau.alpha2=1,y.T=y.Timp)


# Build the model.
TD_model <- nimbleModel(model_code, TD_constants,TD_data,TD_inits)
TD_compiled_model <- compileNimble(TD_model,resetFunctions = TRUE)

# Set up samplers.
TD_mcmc_conf <- configureMCMC(TD_model,monitors=c("beta0.T","beta0.D","beta","tau.V1","tau.V2","alpha","alpha2","tau.alpha",
                                                  "U","aa","aa2","V1","V2"),useConjugacy = TRUE)

TD_mcmc<-buildMCMC(TD_mcmc_conf)
TD_compiled_mcmc<-compileNimble(TD_mcmc, project = TD_model,resetFunctions = TRUE)

# Run the model 
MM = 100000
st<-Sys.time()
TD_samples=runMCMC(TD_compiled_mcmc,inits=TD_inits,
                   nchains = 1, nburnin=MM/2,niter = MM,samplesAsCodaMCMC = TRUE,thin=5,
                   summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st

#save(TD_samples,file="FactorModelOutput.Rda")


samples <- as.data.frame(TD_samples)
U_shared = rep(0,n)
for(i in 1:88){
  U_shared[i]= sum(samples[i])/10000
}

Loadings = rep(0,n)
for(i in 1:88){
  Loadings[i]= sum(samples[i+266])/10000
}

Loadings2 = rep(0,n)
for(i in 1:88){
  Loadings2[i]= sum(samples[i+354])/10000
}

V1 = rep(0,n)
for(i in 1:88){
  V1[i]= sum(samples[i+88])/10000
}

V2 = rep(0,n)
for(i in 1:88){
  V2[i]= sum(samples[i+176])/10000
}

summary(samples[,354:355])
sd(samples$`beta[1]`)
sd(samples$`beta[2]`)

