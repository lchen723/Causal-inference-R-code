#####################################
#######  Data Generation  ###########
#####################################

library(ncvreg)
library(MASS)
library(glmnet)
#### General Settings
GAMMA = NULL
DATA = NULL

n = 400       # sample size    ## warnings would happen if n is "small" enough, e.g., n = 100
m = 500        # num of repitition of simulation
px = 15      # dim of parameters
pz = 15
p = px + pz

#### process of parameter
gamma_0 = 1


gamma_X = c(rep(1,ceiling(px/4)),rep(-1,ceiling(px/4)) , rep(0,px-2*ceiling(px/4)))
gamma_Z = c(rep(1,ceiling(pz/4)),rep(-1,ceiling(pz/4)) , rep(0,pz-2*ceiling(pz/4)))


gamma = c(gamma_X, gamma_Z)

#### process of covariate

mu_X = rep(0,px)
mu_Z = rep(0,pz)

Sigma_X = matrix(0,px,px)
Sigma_Z = matrix(0,pz,pz)

sigma_X = 1
sigma_Z = 1
Sigma_e = diag(0.2,px)
rho_X = 0.5
rho_Z = 0.5

for(i in 1:px)
{
for(j in 1:pz)
{
Sigma_X[i,j] = sigma_X * rho_X^(abs(i-j))
Sigma_Z[i,j] = sigma_Z * rho_Z^(abs(i-j))
}
}

for(time in 1:m) # repitition starts
{

X = mvrnorm(n, mu_X, Sigma_X, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

Z = mvrnorm(n, mu_Z, Sigma_Z, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

#### measurement error
err = mvrnorm(n, mu_X, Sigma_e , tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
W = X + err

GX = rep(gamma_X,n)
GZ = rep(gamma_Z,n)
#t(matrix(GX,px,n)) * X  # X^T gamma_X
#t(matrix(GZ,pz,n)) * Z  # Z^T gamma_Z

g_link = X%*%gamma_X + Z%*%gamma_Z
e = 1/(1+exp(-g_link))         # pass through an inv-logit function
#e =  1 - exp(-exp(g_link))     # complement log-log model form
T = rbinom(n,1,e)
  
#### M1
Y = T + abs(X%*%gamma_X + Z%*%gamma_Z) + rnorm(n,0,1)
#######

#### M2
#T = rbinom(n,1,e)
#pi = 1/(1+exp(-(T+g_link)))
#Y = rbinom(n,1,pi)
#######  
  
  


#data = cbind(T,X,Z)    # true data frame
data = cbind(Y,T,W,Z)    # naive data frame

DATA[[time]] = data



}

####    End of Data Generarion   ####


