#############################
####   SIMEX Algorithm   ####
#############################


######  Function for SIMEX Setup

library(glmnet)
library(MASS)

Psi = seq(0,2,length=100)
K = seq(1,100,1)


SIMEX_S1 = function(data, psi, k) {   # Step 1 in the SIMEX algorithm
 
## Algorithm settings

 Y = data[,1]
 T = data[,2]             
 W = data[,3:(px+2)]
 Z = data[,(px+3):(p+2)]
#set.seed(k)
 e = mvrnorm(n,  mu_X, Sigma_e, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
 W_S = W + sqrt(psi) * e

## End of settings

 x = cbind(W_S,Z)
 t = T
 Da = cbind(t,x)
 Da = data.frame(Da)
 output = glm(t~., data = Da,family = "binomial")$coef[-1]

 return(output)
}

SIMEX_S2 = function(data,Psi,K) {   # Step 2 in SIMEX algorithm
GAMMA_Psi = NULL

for(psi in 1 : length(Psi)) {
G_coll = NULL                   # collection for all k in the following loop

for(k in 1 : length(K)) {

GAMMA_psi_k = SIMEX_S1(data,Psi[psi],K[k])
G_coll = rbind(G_coll, GAMMA_psi_k)

}
GAMMA_Psi[[psi]] = G_coll

}

return(GAMMA_Psi)

}






SIMEX_S3 = function(data, Psi, K) {  # Step 3 in SIMEX algorithm

result = SIMEX_S2(data,Psi,K)

vt = NULL
for(psi in 1:length(Psi)) {

vt = rbind(vt,colMeans(result[[psi]]))

}

B = cbind(1,Psi,Psi^2)

beta = solve(t(B)%*%B) %*% (t(B)%*%vt)

ext = c(1,-1,1)

output = ext %*%beta

return(output)

}

#### End of presenting functions 


####  Simulation for proposed method

GAMMA_SIMEX = NULL

for(time in 1:m)
{

data = DATA[[time]]

test = SIMEX_S3(data,Psi,K)

GAMMA_SIMEX = rbind(GAMMA_SIMEX,test)

}



####################################

y = colMeans(GAMMA_SIMEX)
V = diag(1,length(y),length(y))    # Yi mentioned that diagonal matrix is one of the candidate

##### Step 4 with Lasso Approach

LASSO = function(V,y,method="lasso") {   # y is SIMEX estimator and V is weighted matrix

#### glmnet package
#var_sel = glmnet(V, y , family = "gaussian", alpha = 1)
#lambda = which.max(obj2)       # search max
#lambda = which.max(obj2[7:length(obj2)])+6   # search max
#coef = coef[, lambda]
#obj2 = -dev2 + log(n) * reg.df2   # - 2*logL + 2*log(n) x df  in Step 4 by Yi and Chen (201x)
#########


#### ncvreg package   # for lasso, both results are similar.
if(method =="lasso"){
var_sel = ncvreg(V, y , family = "gaussian",penalty="lasso" ,alpha = 1) }

if(method =="scad"){
var_sel = ncvreg(V, y , family = "gaussian",penalty="SCAD" ,alpha = 1)}

coef = rbind(var_sel$a0,as.matrix(var_sel$beta))  # extract coefficients at all values of lambda,  including the intercept
dev2 = deviance(var_sel)
reg.df2 = var_sel$df


obj2 = BIC(var_sel)

lambda = which.min(obj2)       # search max
coef = coef[ , lambda]




output = coef[-1]
return(output)


}



est_lasso = LASSO(V,y,method="lasso")
est_scad = LASSO(V,y,method="scad")




