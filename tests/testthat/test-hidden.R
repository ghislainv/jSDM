#context("test-hidden")

# Data simulation
#= Number of sites
nsite <- 50
#= Number of species
nsp <- 5
#= Set seed for repeatability
seed <- 1234
#= Number of visits associated to each site
set.seed(seed)
visits <- rpois(nsite,3)
visits[visits==0] <- 1
#= Number of latent variables
n_latent <- 2

form.Tr <- function(trait_formula, trait_data,X){
  data <- trait_data
  # add column of 1 with names of covariables in site_data 
  data[,colnames(X)] <- 1
  mf.suit.tr <- model.frame(formula=trait_formula, data=data)
  # full design matrix corresponding to formula
  mod.mat <- model.matrix(attr(mf.suit.tr,"terms"), data=mf.suit.tr)
  # Remove duplicated columns to get design matrix for traits 
  Tr <- as.matrix(mod.mat[,!duplicated(mod.mat,MARGIN=2)])
  colnames(Tr) <- colnames(mod.mat)[!duplicated(mod.mat,MARGIN=2)]
  # Rename columns according to considered trait 
  for(p in 1:np){
    if(sum(colnames(Tr)==colnames(X)[p])==0){
      colnames(Tr) <- gsub(pattern=paste0(":",colnames(X)[p]), replacement="",
                           x=colnames(Tr), fixed=TRUE)
      colnames(Tr) <- gsub(pattern=paste0(colnames(X)[p],":"), replacement="",
                           x=colnames(Tr), fixed=TRUE)
    }
  }
  nt <- ncol(Tr)
  n_Tint <- sum(sapply(apply(Tr,2,unique), FUN=function(x){all(x==1)}))
  col_Tint <- which(sapply(apply(Tr,2,unique), FUN=function(x){all(x==1)}))
  gamma_zeros <- matrix(0,nt,np)
  rownames(gamma_zeros) <- colnames(Tr)
  colnames(gamma_zeros) <- colnames(X)
  for(t in 1:nt){
    for(p in 1:np){
      term <- c(grep(paste0(colnames(X)[p],":"), colnames(mod.mat), value=TRUE, fixed=TRUE),grep(paste0(":",colnames(X)[p]), colnames(mod.mat), value=TRUE, fixed=TRUE))
      if(length(term)==0) next
      # fixed=TRUE pattern is a string to be matched as is 
      # not a regular expression because of special characters in formula (^, /, [, ...) 
      gamma_zeros[t,p] <- length(c(grep(paste0(":",colnames(Tr)[t]), term, fixed=TRUE),grep(paste0(colnames(Tr)[t],":"), term, fixed=TRUE)))
    }
    gamma_zeros[t,1] <- length(which(colnames(mod.mat)==colnames(Tr)[t]))  
  }
  gamma_zeros[col_Tint,] <- 1
  return(list(gamma_zeros=gamma_zeros,Tr=Tr))
}

# Ecological process (suitability)
x1 <- rnorm(nsite,0,1)
set.seed(2*seed)
x2 <- rnorm(nsite,0,1)
X <- cbind(rep(1,nsite),x1,x2)
colnames(X) <- c("Int", "x1", "x2")
np <- ncol(X)
trait_data <- data.frame(WSD=scale(runif(nsp,0,1000)), SLA=scale(runif(nsp,0,250)))
trait_formula <- ~ WSD + SLA + x1:I(WSD^2) + I(x1^2):SLA + x2:I(SLA^2) + I(x2^2):WSD
result <- form.Tr(trait_formula,trait_data,X)
Tr <- result$Tr
nt <- ncol(Tr)
gamma_zeros <- result$gamma_zeros
gamma.target <- matrix(runif(nt*np,-2,2), byrow=TRUE, nrow=nt)
#= Species effect beta
mu_beta <- as.matrix(Tr) %*% (gamma.target*gamma_zeros)
V_beta <- diag(1,np)
beta.target <- matrix(NA,nrow=np,ncol=nsp)
for(j in 1:nsp){
  beta.target[,j] <- MASS::mvrnorm(n=1, mu=mu_beta[j,], Sigma=V_beta)
}
#= Latent variables W
W <- matrix(rnorm(nsite*n_latent,0,1), nsite, n_latent)
#= Factor loading lambda
lambda.target <- matrix(0, n_latent, nsp)
mat <- t(matrix(runif(nsp*n_latent, -2, 2), byrow=TRUE, nrow=nsp))
lambda.target[upper.tri(mat, diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
#= site effect alpha 
Valpha.target <- 0.5
alpha.target <- rnorm(nsite,0,sqrt(Valpha.target))

#### Poisson ############
log.theta <- X %*% beta.target + W%*%lambda.target + alpha.target
theta <- exp(log.theta)
set.seed(seed)
Y.pois <- apply(theta, 2, rpois, n=nsite) 

#### Binomial ####
#= Number of visits associated to each site
visits <- rpois(nsite,3)
visits[visits==0] <- 1
logit.theta <- X %*% beta.target + W%*%lambda.target + alpha.target
theta <- jSDM::inv_logit(logit.theta)
set.seed(seed)
Y.bin <- apply(theta, 2, rbinom, n=nsite, size=visits)

#### Binomial long format 

# Ecological process (suitability)
## X
x1 <- rnorm(nsite,0,1)
x1.2 <- scale(x1^2)
X <- cbind(rep(1,nsite),x1,x1.2)
colnames(X) <- c("Int","x1","x1.2")
np <- ncol(X)
## W
W <- matrix(rnorm(nsite*n_latent,0,1),nrow=nsite,byrow=TRUE)
## D
SLA <- runif(nsp,-1,1)
D <- data.frame(x1.SLA= scale(c(x1 %*% t(SLA))))
nd <- ncol(D)
## parameters
beta.target <- t(matrix(runif(nsp*np,-2,2), byrow=TRUE, nrow=nsp))
mat <- t(matrix(runif(nsp*n_latent,-2,2), byrow=TRUE, nrow=nsp))
diag(mat) <- runif(n_latent,0,2)
lambda.target <- matrix(0,n_latent,nsp)
gamma.target <-runif(nd,-1,1)
# constraints of identifiability
beta.target[,1] <- 0.0
lambda.target[upper.tri(mat,diag=TRUE)] <- mat[upper.tri(mat, diag=TRUE)]
#= Variance of random site effect
V_alpha.target <- 0.5
#= Random site effect
alpha.target <- rnorm(nsite,0,sqrt(V_alpha.target))
## probit_theta
probit_theta <- c(X %*% beta.target) + c(W %*% lambda.target) + as.matrix(D) %*% gamma.target + rep(alpha.target,nsp)
# Supplementary observation (each site have been visited twice)
# Environmental variables at the time of the second visit
x1_supObs <- rnorm(nsite,0,1)
x1.2_supObs <- scale(x1^2)
X_supObs <- cbind(rep(1,nsite),x1_supObs,x1.2_supObs)
D_supObs <- data.frame(x1.SLA=scale(c(x1_supObs %*% t(SLA))))
probit_theta_supObs <- c(X_supObs%*%beta.target) + c(W%*%lambda.target) + as.matrix(D_supObs) %*% gamma.target + alpha.target
probit_theta <- c(probit_theta, probit_theta_supObs)
nobs <- length(probit_theta)
e <- rnorm(nobs,0,1)
Z_true <- probit_theta + e
Y<-rep(0,nobs)
for (n in 1:nobs){
  if ( Z_true[n] > 0) {Y[n] <- 1}
}
Id_site <- rep(1:nsite,nsp)
Id_sp <- rep(1:nsp,each=nsite)
data <- data.frame(site=rep(Id_site,2), species=rep(Id_sp,2), Y=Y,
                   x1=c(rep(x1,nsp),rep(x1_supObs,nsp)), x1.2=c(rep(x1.2,nsp),rep(x1.2_supObs,nsp)),
                   x1.SLA=c(D$x1.SLA,D_supObs$x1.SLA))
# missing observation
data <- data[-1,]
nobs <- nobs-1


lambda_start <- matrix(1, n_latent, nsp)
lambda_start[lower.tri(lambda_start)] <- 0
test_that("checks works", {
  expect_equal(check.mcmc.parameters(1000, 1000, 1), 0)
  expect_equal(check.verbose(1), 0)
  expect_equal(check.verbose(0), 0)
  expect_equal(check.Y.poisson(Y.pois), 0)
  expect_equal(check.Y.binomial(Y.bin, T=visits), 0)
  expect_equal(check.T.binomial(visits, nsite), 0)
  expect_equal(check.X(X, nsite), 0)
  expect_equal(form.beta.start(0,np), rep(0,np))
  expect_equal(form.beta.start(rep(0,np),np), rep(0,np))
  expect_equal(form.beta.start.sp(0, np, nsp), matrix(0, np, nsp))
  expect_equal(form.lambda.start.sp(0, n_latent, nsp), diag(1, n_latent, nsp))
  expect_equal(form.lambda.start.sp(1, n_latent, nsp), lambda_start)
  expect_equal(form.W.start.sp(0, nsite, n_latent), matrix(0, nsite, n_latent))
  expect_equal(form.alpha.start.sp(0, nsite), rep(0, nsite))
  expect_equal(form.b.start(0, nd), rep(0, nd))
  expect_equal(form.gamma.start.mat(0, np, nt), matrix(0, np, nt))
  expect_equal(check.mub(0, nd), rep(0, nd))
  expect_equal(check.Vb.mat(1, nd), diag(rep(1,nd)))
  expect_equal(check.Vb.mat(matrix(1,nd,nd), nd), matrix(1,nd, nd))
  expect_equal(check.mugamma.mat(0, nt, np), matrix(0, nt, np))
  expect_equal(check.Vgamma.mat(1, nt, np), matrix(1, nt, np))
  expect_equal(check.Vgamma.mat(1:nt, nt, np), matrix(1:nt, nt, np))
  expect_equal(check.Vgamma.mat(1:np, nt, np), matrix(1:np, nt, np, byrow=TRUE))
  expect_equal(check.Vgamma.mat(matrix(1:(np*nt), nt, np), nt, np), matrix(1:(np*nt), nt, np))
  expect_equal(check.Valpha(1), 1)
  expect_equal(check.Vlambda(1, n_latent), rep(1,n_latent))
  expect_equal(check.Vlambda(rep(1,n_latent), n_latent), rep(1,n_latent))
  expect_equal(check.mulambda(0, n_latent), rep(0, n_latent))
  expect_equal(check.mulambda(rep(0,n_latent), n_latent), rep(0,n_latent))
  expect_equal(check.Vbeta(1, np), rep(1,np))
  expect_equal(check.Vbeta(rep(1,np), np), rep(1,np))
  expect_equal(check.mubeta(0, np), rep(0, np))
  expect_equal(check.mubeta(rep(0,np), np), rep(0, np))
  expect_equal(check.Vbeta.mat(1,np), diag(rep(1,np)))
  expect_equal(check.Vbeta.mat(rep(1,np),np), diag(rep(1,np)))
  expect_equal(check.Vbeta.mat(diag(rep(1,np)),np), diag(rep(1,np)))
  expect_equal(check.Vlambda.mat(1, n_latent), diag(rep(1,n_latent)))
  expect_equal(check.Vlambda.mat(rep(1, n_latent), n_latent), diag(rep(1, n_latent)))
  expect_equal(check.Vlambda.mat(diag(rep(1, n_latent)), n_latent), diag(rep(1, n_latent)))
})

