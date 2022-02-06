#################################################
####   SIMEX Algorithm: estimation for ATE   ####
#################################################


tau_S1 = function(data, psi, k, gamma) {   # Step 1 in the SIMEX algorithm
 
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
 g_link = x%*%gamma
 ps = 1/(1+exp(-g_link))
 output = (sum(T*Y/ps) / sum(T/ps)) - (sum((1-T)*Y/(1-ps)) / sum((1-T)/(1-ps)) )

 return(output)
}

tau_S2 = function(data,Psi,K, gamma) {   # Step 2 in SIMEX algorithm
GAMMA_Psi = NULL

for(psi in 1 : length(Psi)) {
G_coll = NULL                   # collection for all k in the following loop

for(k in 1 : length(K)) {

GAMMA_psi_k = tau_S1(data,Psi[psi],K[k], gamma)
G_coll = rbind(G_coll, GAMMA_psi_k)

}
GAMMA_Psi[[psi]] = G_coll

}

return(GAMMA_Psi)

}


#test = SIMEX_S2(data,Psi,K, gamma)



tau_S3 = function(data, Psi, K, gamma) {  # Step 3 in SIMEX algorithm

result = tau_S2(data,Psi,K, gamma)

vt = NULL
for(psi in 1:length(Psi)) {

vt = rbind(vt,colMeans(result[[psi]]))

}

 
##  quadratic extrapolant
B = cbind(1,Psi,Psi^2)
beta = solve(t(B)%*%B) %*% (t(B)%*%vt)
ext = c(1,-1,1)

##  linear extrapolant
#B = cbind(1,Psi)
#beta = solve(t(B)%*%B) %*% (t(B)%*%vt)
#ext = c(1,-1)
 
 
output = ext %*%beta

return(output)

}

tau_lasso = NULL
tau_scad = NULL

for(timetau in 1:time)
{

data = DATA[[timetau]]

tauL = tau_S3(data, Psi, K, gamma = est_lasso)
tauS = tau_S3(data, Psi, K, gamma = est_scad)

tau_lasso = rbind(tau_lasso,tauL)
tau_scad = rbind(tau_scad,tauS)

}

mean(tau_lasso)
mean(tau_scad)


