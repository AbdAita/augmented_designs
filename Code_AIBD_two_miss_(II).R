mu.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,treatment))
  m=ncol(x2)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  for (l in 1:b) {  
    SBT[l]<-t(x1[,l])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:m) {  
    TCA[j]<-sum(x2[,j]*x4[,1]*y)
  }
  rr<-data.frame(ee, TCA)
  r <- with(rr, which(TCA==0, arr.ind=TRUE))
  u<-x4[,2]*y
  l<-data.frame(obs, u)
  ll <- with(l, which(u==0, arr.ind=TRUE))
  newd1 <- l[-ll, ]
  newd <- rr[-r, ]
  f<-x4[,2]*x1
  h<-colSums(f)
  yBM..<-newd[,2]
  yBMjg<-newd1[,2]
  c=length(yBM..)
  k1=colSums(x1)
  N=sum(k)
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  yBi..=BC
  y.M..=t(x2)%*%y
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  i=n_block_i
  ii=n_block_ii
  j=n_check_j
  jj=n_check_jj
  #i=1
  #ii=3
  #j=3
  #jj=3
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
  beta=c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  mu=(y....-((b-1)*(sum(yBM..)-yBM..[j])/b)-((b-3)*(yBM..[j])/(b-2))-sum(nn_new*beta)+nn_new[i]*beta[i]+nn_new[ii]*beta[ii]+((1-(b-2)*nn_new[i])*beta[i])/(b-2)+((1-(b-2)*nn_new[ii])*beta[ii])/(b-2))/(n_t+c)
  mu_RT=(y....-((nn_new[i]-1)*yB...[i])/(c+nn_new[i]-1)-((nn_new[ii]-1)*yB...[ii])/(c+nn_new[ii]-1)-sum((nn_new*yB...)/(c+nn_new))+((nn_new[i])*yB...[i])/(c+nn_new[i])+((nn_new[ii])*yB...[ii])/(c+nn_new[ii]))/(b*c)
  mu_RB=(y....-((b-1)*(sum(yBM..)-yBM..[j]))/b-((b-3)*(yBM..[j])/(b-2)))/(n_t+c)
  anova.table=cbind(mu_q,mu,mu_RT,mu_RB)
  colnames(anova.table)=c("mu_q","mu","mu_RT","mu_RB")
  return(anova.table)
} 
################################################################
Beta.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,treatment))
  m=ncol(x2)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  for (l in 1:b) {  
    SBT[l]<-t(x1[,l])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:m) {  
    TCA[j]<-sum(x2[,j]*x4[,1]*y)
  }
  rr<-data.frame(ee, TCA)
  r <- with(rr, which(TCA==0, arr.ind=TRUE))
  u<-x4[,2]*y
  l<-data.frame(obs, u)
  ll <- with(l, which(u==0, arr.ind=TRUE))
  newd1 <- l[-ll, ]
  newd <- rr[-r, ]
  f<-x4[,2]*x1
  h<-colSums(f)
  yBM..<-newd[,2]
  yBMjg<-newd1[,2]
  c=length(yBM..)
  k1=colSums(x1)
  N=sum(k)
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  yBi..=BC
  y.M..=t(x2)%*%y
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  i=n_block_i
  ii=n_block_ii
  j=n_check_j
  jj=n_check_jj
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
  beta=c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  mu_RT=(y....-((nn_new[i]-1)*yB...[i])/(c+nn_new[i]-1)-((nn_new[ii]-1)*yB...[ii])/(c+nn_new[ii]-1)-sum((nn_new*yB...)/(c+nn_new))+((nn_new[i])*yB...[i])/(c+nn_new[i])+((nn_new[ii])*yB...[ii])/(c+nn_new[ii]))/(b*c)
  beta_RT=c()
  for (a in 1:b) { 
    if(a==i){
      beta_RT[a]<-(yB...[i])/(c+nn_new[a]-1)-mu_RT
    }
    if(a==ii){
      beta_RT[a]<-(yB...[ii])/(c+nn_new[a]-1)-mu_RT
    }
    if(a!=i & a!=ii){
      beta_RT[a]<-(yB...[a])/(c+nn_new[a])-mu_RT
    }
  }
  anova.table=cbind(beta,beta_RT)
  colnames(anova.table)=c("beta","beta_RT")
  return(anova.table)
} 
################################################################
Check.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,treatment))
  m=ncol(x2)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  for (l in 1:b) {  
    SBT[l]<-t(x1[,l])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:m) {  
    TCA[j]<-sum(x2[,j]*x4[,1]*y)
  }
  rr<-data.frame(ee, TCA)
  r <- with(rr, which(TCA==0, arr.ind=TRUE))
  u<-x4[,2]*y
  l<-data.frame(obs, u)
  ll <- with(l, which(u==0, arr.ind=TRUE))
  newd1 <- l[-ll, ]
  newd <- rr[-r, ]
  f<-x4[,2]*x1
  h<-colSums(f)
  yBM..<-newd[,2]
  yBMjg<-newd1[,2]
  c=length(yBM..)
  k1=colSums(x1)
  N=sum(k)
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  yBi..=BC
  y.M..=t(x2)%*%y
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  i=n_block_i
  ii=n_block_ii
  j=n_check_j
  jj=n_check_jj
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
  beta=c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  mu=(y....-((b-1)*(sum(yBM..)-yBM..[j])/b)-((b-3)*(yBM..[j])/(b-2))-sum(nn_new*beta)+nn_new[i]*beta[i]+nn_new[ii]*beta[ii]+((1-(b-2)*nn_new[i])*beta[i])/(b-2)+((1-(b-2)*nn_new[ii])*beta[ii])/(b-2))/(n_t+c)
  taucheck<-c()
  for (e in 1:c) { 
    if(e==j){
      taucheck[e]<-(yBM..[j]+beta[i]+beta[ii])/(b-2)-mu
    }
    if(e!=j){
      taucheck[e]<-(yBM..[e])/b-mu
    }
  }
  mu_RB=(y....-((b-1)*(sum(yBM..)-yBM..[j]))/b-((b-3)*(yBM..[j])/(b-2)))/(n_t+c)
  taucheck_NB<-c()
  for (e in 1:c) { 
    if(e==j){
      taucheck_NB[e]<-yBM..[j]/(b-2)-mu_RB
    }
    if(e!=j){
      taucheck_NB[e]<-yBM..[e]/b-mu_RB
    }
  } 
  xx22=as.matrix(table(obs,1:N))
  CC<-t(x2)%*%x4[,1]
  CCC=CC[apply(CC!=0, 1, all),]
  check_tretment=CCC
  Rep_Check=cbind(check_tretment)
  anova.table=data.frame(Rep_Check,taucheck,taucheck_NB)
  colnames(anova.table)=c("Rep_Check","taucheck","taucheck_NB")
  return(anova.table)
} 
################################################################
New_treatment.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,treatment))
  m=ncol(x2)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  for (l in 1:b) {  
    SBT[l]<-t(x1[,l])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:m) {  
    TCA[j]<-sum(x2[,j]*x4[,1]*y)
  }
  rr<-data.frame(ee, TCA)
  r <- with(rr, which(TCA==0, arr.ind=TRUE))
  u<-x4[,2]*y
  l<-data.frame(obs, u)
  ll <- with(l, which(u==0, arr.ind=TRUE))
  newd1 <- l[-ll, ]
  newd <- rr[-r, ]
  f<-x4[,2]*x1
  h<-colSums(f)
  yBM..<-newd[,2]
  yBMjg<-newd1[,2]
  c=length(yBM..)
  k1=colSums(x1)
  N=sum(k)
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  yBi..=BC
  y.M..=t(x2)%*%y
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  i=n_block_i
  ii=n_block_ii
  j=n_check_j
  jj=n_check_jj
  nn_new=t(x1)%*%x4[,2]
  nn_new[1:2]
  n_new=nn_new[1]
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
  beta=c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  mu=(y....-((b-1)*(sum(yBM..)-yBM..[j])/b)-((b-3)*(yBM..[j])/(b-2))-sum(nn_new*beta)+nn_new[i]*beta[i]+nn_new[ii]*beta[ii]+((1-(b-2)*nn_new[i])*beta[i])/(b-2)+((1-(b-2)*nn_new[ii])*beta[ii])/(b-2))/(n_t+c)
  xx22=as.matrix(table(obs,1:N))
  BB<-t(xx22)%*%y
  AA<-t(xx22)%*%x4[,2]
  M<-AA*BB
  MM=M[apply(M!=0, 1, all),] 
  ff<-c()
  nnn<-0
  for(iii in 1:b){
    for(jjj in nnn+1:nn_new[iii]){ 
      ff[jjj]=MM[jjj]-beta[iii]-mu
    }
    nnn<-nnn+nn_new[iii]
  }
  ni<-0
  test_tretment<-c()
  for(iii in 1:N){
    if(AA[iii]==1){
      ni<-ni+1
      test_tretment[ni]<-treatment[iii]
    }
  }
  New_Treatment=cbind(test_tretment)
  taunew=cbind(ff)
  mu_RB=(y....-((b-1)*(sum(yBM..)-yBM..[j]))/b-((b-3)*(yBM..[j])/(b-2)))/(n_t+c)
  tau_new_RB<-c()
  for (h in 1:n_t) {
    tau_new_RB[h]=MM[h]-mu_RB
  }
  taunew_RB=cbind(tau_new_RB)
  anova.table=data.frame(New_Treatment,taunew,taunew_RB)
  colnames(anova.table)=c("New_Treatment","taunew","taunew_RB")
  return(anova.table)
} 
################################################################
aov1.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,treatment))
  m=ncol(x2)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  for (l in 1:b) {  
    SBT[l]<-t(x1[,l])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:m) {  
    TCA[j]<-sum(x2[,j]*x4[,1]*y)
  }
  rr<-data.frame(ee, TCA)
  r <- with(rr, which(TCA==0, arr.ind=TRUE))
  u<-x4[,2]*y
  l<-data.frame(obs, u)
  ll <- with(l, which(u==0, arr.ind=TRUE))
  newd1 <- l[-ll, ]
  newd <- rr[-r, ]
  f<-x4[,2]*x1
  h<-colSums(f)
  yBM..<-newd[,2]
  yBMjg<-newd1[,2]
  c=length(yBM..)
  k1=colSums(x1)
  N=sum(k)
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  yBi..=BC
  y.M..=t(x2)%*%y
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  i=n_block_i
  ii=n_block_ii
  j=n_check_j
  jj=n_check_jj
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
  beta=c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  mu=(y....-((b-1)*(sum(yBM..)-yBM..[j])/b)-((b-3)*(yBM..[j])/(b-2))-sum(nn_new*beta)+nn_new[i]*beta[i]+nn_new[ii]*beta[ii]+((1-(b-2)*nn_new[i])*beta[i])/(b-2)+((1-(b-2)*nn_new[ii])*beta[ii])/(b-2))/(n_t+c)
  mu_RT=(y....-((nn_new[i]-1)*yB...[i])/(c+nn_new[i]-1)-((nn_new[ii]-1)*yB...[ii])/(c+nn_new[ii]-1)-sum((nn_new*yB...)/(c+nn_new))+((nn_new[i])*yB...[i])/(c+nn_new[i])+((nn_new[ii])*yB...[ii])/(c+nn_new[ii]))/(b*c)
  mu_RB=(y....-((b-1)*(sum(yBM..)-yBM..[j]))/b-((b-3)*(yBM..[j])/(b-2)))/(n_t+c)
  beta=c()
  beta_RT=c()
  taucheck<-c()
  taucheck_NB<-c()
  test_tretment<-c()
  tau_new_RB<-c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  for (a in 1:b) { 
    if(a==i){
      beta_RT[a]<-(yB...[i])/(c+nn_new[a]-1)-mu_RT
    }
    if(a==ii){
      beta_RT[a]<-(yB...[ii])/(c+nn_new[a]-1)-mu_RT
    }
    if(a!=i & a!=ii){
      beta_RT[a]<-(yB...[a])/(c+nn_new[a])-mu_RT
    }
  }
  
  for (e in 1:c) { 
    if(e==j){
      taucheck[e]<-(yBM..[j]+beta[i]+beta[ii])/(b-2)-mu
    }
    if(e!=j){
      taucheck[e]<-(yBM..[e])/b-mu
    }
  }
  for (e in 1:c) { 
    if(e==j){
      taucheck_NB[e]<-yBM..[j]/(b-2)-mu_RB
    }
    if(e!=j){
      taucheck_NB[e]<-yBM..[e]/b-mu_RB
    }
  }
  
  xx22=as.matrix(table(obs,1:N))
  BB<-t(xx22)%*%y
  AA<-t(xx22)%*%x4[,2]
  M<-AA*BB
  MM=M[apply(M!=0, 1, all),] 
  ff<-c()
  nnn<-0
  for(iii in 1:b){
    for(jjj in nnn+1:nn_new[iii]){ 
      ff[jjj]=MM[jjj]-beta[iii]-mu
    }
    nnn<-nnn+nn_new[iii]
  }
  ni<-0
  for(iii in 1:N){
    if(AA[iii]==1){
      ni<-ni+1
      test_tretment[ni]<-treatment[iii]
    }
  }
  New_Treatment=cbind(test_tretment)
  n_treatment=ff
  taunew=cbind(ff)
  for (h in 1:n_t) {
    tau_new_RB[h]=MM[h]-mu_RB
  }
  taunew_RB=cbind(tau_new_RB)
  #################################
  R_full=y....*mu+sum(yB...*beta)+sum(yBM..*taucheck)+sum(yBMjg*n_treatment)
  R_reduced_B=y....*mu_RB+sum(yBM..*taucheck_NB)+sum(yBMjg*tau_new_RB)
  R_reduced_Tr=y....*mu_RT+sum(yB...*beta_RT)
  #################################
  blockSSadj=R_full- R_reduced_B
  treatmentSSadj=R_full- R_reduced_Tr
  blockSSunadj=sum(yB...^2/k1)-y....^2/(N)
  totalSS=t(y)%*%y-y....^2/(N)
  errorSS=totalSS- blockSSunadj-treatmentSSadj
  #############################
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
    tauCheck<-c()
    for (e in 1:c) { 
      if(e==j){
        tauCheck[e]<-(yBM..[j]+beta[i]+beta[ii])/(b-2)-mu_q
      }
      if(e!=j){
        tauCheck[e]<-(yBM..[e])/b-mu_q
      }
    }
  mq_RT=(sumcheck+(yBi..[i]++yBi..[ii])/(c-1))/(c*b)
  betaCheck_RT<-c()
  for (a in 1:b) { 
    if(a==i){
      betaCheck_RT[a]<-yBi..[i]/(c-1)-mq_RT
    }
    if(a==ii){
      betaCheck_RT[a]<-yBi..[ii]/(c-1)-mq_RT
    }
    if(a!=i & a!=ii){
      betaCheck_RT[a]<-yBi..[a]/c-mq_RT
    }
  }
  R_full_check=sumcheck*mu_q+sum(yBi..*beta)+sum(yBM..*tauCheck)
  R_reduced_check_Tr=sumcheck*mq_RT+sum(yBi..*betaCheck_RT)
  checkSS=R_full_check- R_reduced_check_Tr
  new_and_new_vs._checkSS=treatmentSSadj-checkSS
  dfE=N-b-m+1 
  meanB=blockSSunadj/(b-1)
  meanT=treatmentSSadj/(m-1)
  meanE=errorSS/dfE
  meanCheck=checkSS/(c-1)
  meanNew_and_new_vs._check=new_and_new_vs._checkSS/(m-c)
  blockF=meanB/meanE
  treatmentF=meanT/meanE
  checkF=meanCheck/meanE
  New_and_new_vs._checkF=meanNew_and_new_vs._check/meanE
  p.value.block=pf(blockF, df1=(b-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.treatment=pf(treatmentF, df1=(m-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.check=pf(checkF, df1=(c-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new_and_new_vs._check=pf(New_and_new_vs._checkF, df1=(m-c), df2=dfE, ncp=0, lower.tail = FALSE)
  df=c(b-1,m-1,c-1,m-c,dfE,N-1)
  ss=c(blockSSunadj,treatmentSSadj,checkSS,new_and_new_vs._checkSS,errorSS,totalSS)
  ms=c(meanB,meanT,meanCheck,meanNew_and_new_vs._check,meanE,NA)
  F=c(NA,treatmentF,checkF,New_and_new_vs._checkF,NA,NA)
  p.value=c(NA,p.value.treatment,p.value.check,p.value.new_and_new_vs._check,NA,NA)
  anova.table=cbind(df,ss,ms,F, p.value)
  rownames(anova.table)=c("Block(Unadj.)","Treatment(Adj.)","  Check(Adj.)","  New & new vs. check","Error","Total")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
} 
###############################################################
###############################################################
aov2.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,block))
  b=ncol(x1)
  x2=as.matrix(table(obs,treatment))
  m=ncol(x2)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  for (l in 1:b) {  
    SBT[l]<-t(x1[,l])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:m) {  
    TCA[j]<-sum(x2[,j]*x4[,1]*y)
  }
  rr<-data.frame(ee, TCA)
  r <- with(rr, which(TCA==0, arr.ind=TRUE))
  u<-x4[,2]*y
  l<-data.frame(obs, u)
  ll <- with(l, which(u==0, arr.ind=TRUE))
  newd1 <- l[-ll, ]
  newd <- rr[-r, ]
  f<-x4[,2]*x1
  h<-colSums(f)
  yBM..<-newd[,2]
  yBMjg<-newd1[,2]
  c=length(yBM..)
  k1=colSums(x1)
  N=sum(k)
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  yBi..=BC
  y.M..=t(x2)%*%y
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  i=n_block_i
  ii=n_block_ii
  j=n_check_j
  jj=n_check_jj
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_q=((b*c-b-2*c)*sumcheck+2*c*(b-2)*yBM..[j]+b*(yBi..[i]+yBi..[ii]))/(b*c*(b-2)*(c-1))
  beta=c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  mu=(y....-((b-1)*(sum(yBM..)-yBM..[j])/b)-((b-3)*(yBM..[j])/(b-2))-sum(nn_new*beta)+nn_new[i]*beta[i]+nn_new[ii]*beta[ii]+((1-(b-2)*nn_new[i])*beta[i])/(b-2)+((1-(b-2)*nn_new[ii])*beta[ii])/(b-2))/(n_t+c)
  mu_RT=(y....-((nn_new[i]-1)*yB...[i])/(c+nn_new[i]-1)-((nn_new[ii]-1)*yB...[ii])/(c+nn_new[ii]-1)-sum((nn_new*yB...)/(c+nn_new))+((nn_new[i])*yB...[i])/(c+nn_new[i])+((nn_new[ii])*yB...[ii])/(c+nn_new[ii]))/(b*c)
  mu_RB=(y....-((b-1)*(sum(yBM..)-yBM..[j]))/b-((b-3)*(yBM..[j])/(b-2)))/(n_t+c)
  beta=c()
  beta_RT=c()
  taucheck<-c()
  taucheck_NB<-c()
  test_tretment<-c()
  tau_new_RB<-c()
  for (a in 1:b) { 
    if(a==i){
      beta[a]<-(b*yBi..[i]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a==ii){
      beta[a]<-(b*yBi..[ii]-sum(yBM..)+yBM..[j])/(b*(c-1))
    }
    if(a!=i & a!=ii){
      beta[a]<-(yBi..[a])/c-mu_q
    }
  }
  for (a in 1:b) { 
    if(a==i){
      beta_RT[a]<-(yB...[i])/(c+nn_new[a]-1)-mu_RT
    }
    if(a==ii){
      beta_RT[a]<-(yB...[ii])/(c+nn_new[a]-1)-mu_RT
    }
    if(a!=i & a!=ii){
      beta_RT[a]<-(yB...[a])/(c+nn_new[a])-mu_RT
    }
  }
  
  for (e in 1:c) { 
    if(e==j){
      taucheck[e]<-(yBM..[j]+beta[i]+beta[ii])/(b-2)-mu
    }
    if(e!=j){
      taucheck[e]<-(yBM..[e])/b-mu
    }
  }
  for (e in 1:c) { 
    if(e==j){
      taucheck_NB[e]<-yBM..[j]/(b-2)-mu_RB
    }
    if(e!=j){
      taucheck_NB[e]<-yBM..[e]/b-mu_RB
    }
  }
  
  xx22=as.matrix(table(obs,1:N))
  BB<-t(xx22)%*%y
  AA<-t(xx22)%*%x4[,2]
  M<-AA*BB
  MM=M[apply(M!=0, 1, all),] 
  ff<-c()
  nnn<-0
  for(iii in 1:b){
    for(jjj in nnn+1:nn_new[iii]){ 
      ff[jjj]=MM[jjj]-beta[iii]-mu
    }
    nnn<-nnn+nn_new[iii]
  }
  ni<-0
  for(iii in 1:N){
    if(AA[iii]==1){
      ni<-ni+1
      test_tretment[ni]<-treatment[iii]
    }
  }
  New_Treatment=cbind(test_tretment)
  n_treatment=ff
  taunew=cbind(ff)
  for (h in 1:n_t) {
    tau_new_RB[h]=MM[h]-mu_RB
  }
  taunew_RB=cbind(tau_new_RB)
  #################################
  R_full=y....*mu+sum(yB...*beta)+sum(yBM..*taucheck)+sum(yBMjg*n_treatment)
  R_reduced_B=y....*mu_RB+sum(yBM..*taucheck_NB)+sum(yBMjg*tau_new_RB)
  R_reduced_Tr=y....*mu_RT+sum(yB...*beta_RT)
  #################################
  blockSSadj=R_full- R_reduced_B
  treatmentSSunadj=sum(y.M..^2/k)-y....^2/(N)
  totalSS=t(y)%*%y-y....^2/(N)
  errorSS=totalSS-blockSSadj-treatmentSSunadj
  checkSS1=(sum(yBM..^2)-yBM..[j]^2)/b+(yBM..[j]^2)/(b-2)-(sumcheck^2/n_c)
  newSS=sum(yBMjg^2)-(sum(yBMjg))^2/n_t
  new_vs._checkSS1=treatmentSSunadj-checkSS1-newSS
  dfE=N-b-m+1 
  meanB=blockSSadj/(b-1)
  meanT=treatmentSSunadj/(m-1)
  meanE=errorSS/dfE
  meanCheck1=checkSS1/(c-1)
  meanNew=newSS/(n_t-1)
  meanNwe_vs._check1=new_vs._checkSS1/1
  blockF=meanB/meanE
  treatmentF=meanT/meanE
  checkF1=meanCheck1/meanE
  newF=meanNew/meanE
  new_vs._checkF1=meanNwe_vs._check1/meanE
  p.value.block=pf(blockF, df1=(b-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.treatment=pf(treatmentF, df1=(m-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.check=pf(checkF1, df1=(c-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new=pf(newF, df1=(n_t-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new_vs._check=pf(new_vs._checkF1, df1=(1), df2=dfE, ncp=0, lower.tail = FALSE)
  df=c(b-1,m-1,c-1,n_t-1,1,dfE,N-1)
  df1=c(b-1,m-1,c-1,n_t-1,1,dfE)
  ss=c(blockSSadj,treatmentSSunadj,checkSS1,newSS,new_vs._checkSS1,errorSS,totalSS)
  ms=c(blockSSadj/(b-1),treatmentSSunadj/(m-1),checkSS1/(c-1),newSS/(n_t-1),new_vs._checkSS1/(1),errorSS/(dfE),NA)
  F=c(blockF,NA,checkF1,newF,new_vs._checkF1,NA,NA)
  p.value=c(p.value.block,NA,p.value.check,p.value.new,p.value.new_vs._check,NA,NA)
  anova.table=cbind(df,ss,ms,F, p.value)
  rownames(anova.table)=c("Block(Adj.)","Treatment(Unadj.)","  Check(Unadj.)","  New","  New vs. check","Error","Total")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
} 
########################################
######################################## 
TMC.II.ARIBD = function(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
{
  output1=mu.ARIBD(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
  output2=Beta.ARIBD(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
  output3=Check.ARIBD(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
  output4=New_treatment.ARIBD(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj) 
  output5=aov1.ARIBD(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
  output6=aov2.ARIBD(obs, block, treatment,type, y,n_block_i,n_block_ii,n_check_j,n_check_jj)
  newlist <- list(output1,output2,output3,output4,output5,output6) 
  return(newlist)
}
#####################################################
