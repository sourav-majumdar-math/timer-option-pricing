rm(list=ls())
gc()

kappa<-2
theta<-0.0324
gamma<-0.1
V_0<-0.0625
rho<- 0.8
V_B<-0.0265
S_0<-100
n_path<-100
Delta<-V_B/n_path
K<-100
r<-0.04

sim_step<-50000

X<-matrix(,nrow = n_path+1,ncol=1)
X_path<-function(kappa,theta,Delta,n_path){
  X[1]=V_0
  for(i in 2:(n_path+1)){
    X[i,1]=X[i-1,1]+kappa*theta*Delta/(X[i-1,1])-kappa*Delta+gamma*rnorm(1,0,sqrt(Delta))
  }
  return(X)
}

tau<-function(n_path,X,Delta)
{
  z<-0
  for(i in 2:(n_path+1)){
    z<-z+Delta/(X[i])
  }
  return(z)
}
a<-function(kappa,theta,gamma,V_0,rho,V_B,tau,V_tau,tau_val){
  return(rho*((V_tau/gamma)-(V_0/gamma))-rho*(kappa*(theta*tau_val-V_B)/gamma)-V_B/2)
}
d_1<-function(r,a_val,rho,V_B,S_0,K,tau_val)
{
  return((log(S_0/K)+r*tau_val+a_val+(1-rho^2)*V_B)/(sqrt((1-rho^2)*V_B)))
}
d_2<-function(d_1_val,rho,V_B){
  return(d_1_val-sqrt((1-rho^2)*V_B))
}

C_mc<-function(S_0,rho,V_B,a_val,d_1_val,K,r,d_2_val,tau_val){
  return(S_0*exp((1-rho^2)*V_B/2)*exp(a_val)*pnorm(d_1_val)-K*exp(-r*tau_val)*pnorm(d_2_val))
}


C<-matrix(,nrow=sim_step,ncol=1)
for(i in 1:sim_step){
  X<-X_path(kappa,theta,Delta,n_path)
  tau_val<-tau(n_path,X,Delta)
  V_tau<-X[(n_path+1)]
  a_val<-a(kappa,theta,gamma,V_0,rho,V_B,tau_val,V_tau,tau_val)
  d_1_val<-d_1(r,a_val,rho,V_B,S_0,K,tau_val)
  d_2_val<-d_2(d_1_val,rho,V_B)
  C[i]<-C_mc(S_0,rho,V_B,a_val,d_1_val,K,r,d_2_val,tau_val)
}

mean(C)
