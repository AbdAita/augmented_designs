aov.AILSD1 = function(obs,row,column,treatment,type, y,n_row,n_column,n_check) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,row))
  b=ncol(x1)
  x2=as.matrix(table(obs,column))
  x3=as.matrix(table(obs,treatment))
  m=ncol(x3)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  SBTT=c()
  for (o in 1:b) {  
    SBT[o]<-t(x1[,o])%*%y
  }
  for (p in 1:b) {  
    SBTT[p]<-t(x2[,p])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  BCC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:b) {  
    BCC[j]<-sum(x2[,j]*x4[,1]*y)
  }
  for (k in 1:m) {  
    TCA[k]<-sum(x3[,k]*x4[,1]*y)
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
  yBi..=BC
  yB.j.=BCC
  yB..k<-newd[,2]
  yBMjg<-newd1[,2]
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  y.M..=t(x2)%*%y
  k1=colSums(x1)
  k2=colSums(x2)
  N=sum(k1)
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  c=length(yB..k)
  i=n_row
  j=n_column
  k=n_check
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_b=((b-3)*sumcheck+yBi..[i]+yB.j.[j]+yB..k[k])/(b*(b-1)*(b-2))
  alpha<-c()
  beta<-c()
  for (a in 1:b) { 
    if(a==i){
  alpha[a]<-((b-1)*yBi..[i]+yB.j.[j]-sum(yB..k)+yB..k[k])/(b*(b-2))
    }
    if(a!=i){
      alpha[a]<-(yBi..[a])/b-mu_b
    }
  }
  for (d in 1:b) { 
    if(d==j){
      beta[d]<-((b-1)*yB.j.[j]+yBi..[i]-sum(yB..k)+yB..k[k])/(b*(b-2))
    }
    if(d!=j){
      beta[d]<-(yB.j.[d])/b-mu_b
    }
  } 
  mu<-(y....-(b-1)*(sum(yB..k)-yB..k[k])/b-((b-2)*yB..k[k])/(b-1)+alpha[i]/(b-1)+beta[j]/(b-1))/(n_t+b)
tau<-c()
  for (e in 1:b) { 
    if(e==k){
      tau[e]<-(yB..k[k]+alpha[i]+beta[j])/(b-1)-mu
    }
    if(e!=k){
      tau[e]<-(yB..k[e])/b-mu
    }
  }  
A<-t(x3[,1:n_t])%*%y
B<-t(x1)%*%x3[,1:n_t]
M=t(B)
C<-t(x2)%*%x3[,1:n_t]
E=t(C)
for (h in 1:n_t) {
  for (v in 1:b) if (M[h,v]!=0 & !is.na(M[h,v]) ) M[h,v]=1
}
for (h in 1:n_t) {
  M[h,]=M[h,]*t(alpha)
}
for (h in 1:n_t) {
  for (v in 1:b) if (E[h,v]!=0 & !is.na(E[h,v]) ) E[h,v]=1
}
for (h in 1:n_t) {
  E[h,]=E[h,]*t(beta)
}
tau_new<-c()
for (h in 1:n_t) {
  tau_new[h]=A[h]-sum(M[h,])-sum(E[h,])-mu
}
R_full=y....*mu+sum(yB...*alpha)+sum(y.M..*beta)+sum(tau_new*A)+sum(tau*yB..k)
alpha_star=((n_t+b^2)*(b+n_new-1)*yB...[i]-(b+n_new-n_t-b^2)*y.M..[j]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
beta_star=((n_t+b^2)*(b+n_new-1)*y.M..[j]-(b+n_new-n_t-b^2)*yB...[i]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
mu_star=(y....+alpha_star+beta_star)/(n_t+b^2-1)
for (a in 1:b) { 
  if(a==i){
    alpha_star[a]<-((n_t+b^2)*(b+n_new-1)*yB...[i]-(b+n_new-n_t-b^2)*y.M..[j]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
  }
  if(a!=i){
    alpha_star[a]<-(yB...[a])/(b+n_new)-mu_star
  }
}
for (d in 1:b) { 
  if(d==j){
    beta_star[d]<-((n_t+b^2)*(b+n_new-1)*y.M..[j]-(b+n_new-n_t-b^2)*yB...[i]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
  }
  if(d!=j){
    beta_star[d]<-(y.M..[d])/(b+n_new)-mu_star
  }
} 
R_reduce=y....*mu_star+sum(yB...*alpha_star)+sum(y.M..*beta_star)
  rowSSunadj=sum(yB...^2/k1)-(y....)^2/N
  columnSSunadj=sum(y.M..^2/k2)-(y....^2/N)
  totalSS=t(y)%*%y-(y....^2/N)
  treatmentSSadj=R_full- R_reduce
  errorSS=totalSS- rowSSunadj-columnSSunadj-treatmentSSadj
  checkSS=sum(yB..k^2)/b+((yBi..[i]+yB.j.[j]+yB..k[k]-sumcheck)^2/(b-2)-(yBi..[i]+yB.j.[j]-sumcheck)^2/(b-1)+(sumcheck*(2*yB..k[k]-sumcheck))/b)/(b-1)
  newSS=sum(yBMjg^2)-(sum(yBMjg))^2/n_t
  new_and_new_vs._checkSS=treatmentSSadj-checkSS
  dfE=N-2*b-m+2 
  meanR=rowSSunadj/(b-1)
  meanC=columnSSunadj/(b-1)
  meanT=treatmentSSadj/(m-1)
  meanE=errorSS/dfE
  meanCheck=checkSS/(c-1)
  meanNew=newSS/(n_t-1)
  meanNew_and_new_vs._check=new_and_new_vs._checkSS/(m-c)
  rowF=meanR/meanE
  columnF=meanC/meanE
  treatmentF=meanT/meanE
  checkF=meanCheck/meanE
  newF=meanNew/meanE
  New_vs._checkSSF=meanNew_vs._checkSS/meanE
  New_and_new_vs._checkF=meanNew_and_new_vs._check/meanE
  p.value.block=pf(rowF, df1=(b-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.block=pf(columnF, df1=(b-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.treatment=pf(treatmentF, df1=(m-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.check=pf(checkF, df1=(c-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new=pf(newF, df1=(n_t-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new_vs._check=pf(New_vs._checkSSF, df1=(m-c), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new_and_new_vs._check=pf(New_and_new_vs._checkF, df1=(m-c), df2=dfE, ncp=0, lower.tail = FALSE)
  df=c(b-1,b-1,m-1,c-1,m-c,dfE,N-1)
  ss=c(rowSSunadj,columnSSunadj,treatmentSSadj,checkSS,new_and_new_vs._checkSS,errorSS,totalSS)
  ms=c(meanR,meanC,meanT,meanCheck,meanNew_and_new_vs._check,meanE,NA)
  F=c(NA,NA,treatmentF,checkF,New_and_new_vs._checkF,NA,NA)
  p.value=c(NA,NA,p.value.treatment,p.value.check,p.value.new_and_new_vs._check,NA,NA)
  anova.table=cbind(df,ss,ms,F, p.value)
  rownames(anova.table)=c("Row(Unadj.)","Column(Unadj.)","Treatment(Adj.)","  Check(Adj.)","  New & new vs. check","Error","Total")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
} 
###############################################################
###############################################################
aov.AILSD2 = function(obs,row,column,treatment,type, y,n_row,n_column,n_check) 
{
  one=matrix(1,length(y),1)
  x1=as.matrix(table(obs,row))
  b=ncol(x1)
  x2=as.matrix(table(obs,column))
  x3=as.matrix(table(obs,treatment))
  m=ncol(x3)
  x4=as.matrix(table(obs,type))
  sumcheck=t(x4[,1])%*%y
  sumtest=t(x4[,2])%*%y
  SBT=c()
  SBTT=c()
  for (o in 1:b) {  
    SBT[o]<-t(x1[,o])%*%y
  }
  for (p in 1:b) {  
    SBTT[p]<-t(x2[,p])%*%y
  }
  ee<-col(x2)[1,]
  k=colSums(x2)
  BC<-c()
  BCC<-c()
  TCA<-c()
  for (i in 1:b) {  
    BC[i]<-sum(x1[,i]*x4[,1]*y)
  }
  for (j in 1:b) {  
    BCC[j]<-sum(x2[,j]*x4[,1]*y)
  }
  for (k in 1:m) {  
    TCA[k]<-sum(x3[,k]*x4[,1]*y)
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
  yBi..=BC
  yB.j.=BCC
  yB..k<-newd[,2]
  yBMjg<-newd1[,2]
  y....=t(one)%*%y
  yB...=t(x1)%*%y
  y.M..=t(x2)%*%y
  k1=colSums(x1)
  k2=colSums(x2)
  S=n_t+1
  SS=n_t+c
  k3=colSums(x3[,S:SS])
  N=sum(k1)
  n_c=colSums(x4)[1]
  n_t=colSums(x4)[2]
  c=length(yB..k)
  i=n_row
  j=n_column
  k=n_check
  nn_new=t(x1)%*%x4[,2]
  n_new=nn_new[1]
  mu_b=((b-3)*sumcheck+yBi..[i]+yB.j.[j]+yB..k[k])/(b*(b-1)*(b-2))
  alpha<-c()
  beta<-c()
  for (a in 1:b) { 
    if(a==i){
      alpha[a]<-((b-1)*yBi..[i]+yB.j.[j]-sum(yB..k)+yB..k[k])/(b*(b-2))
    }
    if(a!=i){
      alpha[a]<-(yBi..[a])/b-mu_b
    }
  }
  for (d in 1:b) { 
    if(d==j){
      beta[d]<-((b-1)*yB.j.[j]+yBi..[i]-sum(yB..k)+yB..k[k])/(b*(b-2))
    }
    if(d!=j){
      beta[d]<-(yB.j.[d])/b-mu_b
    }
  } 
  mu<-(y....-(b-1)*(sum(yB..k)-yB..k[k])/b-((b-2)*yB..k[k])/(b-1)+alpha[i]/(b-1)+beta[j]/(b-1))/(n_t+b)
  tau<-c()
  for (e in 1:b) { 
    if(e==k){
      tau[e]<-(yB..k[k]+alpha[i]+beta[j])/(b-1)-mu
    }
    if(e!=k){
      tau[e]<-(yB..k[e])/b-mu
    }
  }  
  A<-t(x3[,1:n_t])%*%y
  B<-t(x1)%*%x3[,1:n_t]
  M=t(B)
  C<-t(x2)%*%x3[,1:n_t]
  E=t(C)
  for (h in 1:n_t) {
    for (v in 1:b) if (M[h,v]!=0 & !is.na(M[h,v]) ) M[h,v]=1
  }
  for (h in 1:n_t) {
    M[h,]=M[h,]*t(alpha)
  }
  for (h in 1:n_t) {
    for (v in 1:b) if (E[h,v]!=0 & !is.na(E[h,v]) ) E[h,v]=1
  }
  for (h in 1:n_t) {
    E[h,]=E[h,]*t(beta)
  }
  tau_new<-c()
  for (h in 1:n_t) {
    tau_new[h]=A[h]-sum(M[h,])-sum(E[h,])-mu
  }
  R_full=y....*mu+sum(yB...*alpha)+sum(y.M..*beta)+sum(tau_new*A)+sum(tau*yB..k)
  alpha_star=((n_t+b^2)*(b+n_new-1)*yB...[i]-(b+n_new-n_t-b^2)*y.M..[j]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
  beta_star=((n_t+b^2)*(b+n_new-1)*y.M..[j]-(b+n_new-n_t-b^2)*yB...[i]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
  mu_star=(y....+alpha_star+beta_star)/(n_t+b^2-1)
  for (a in 1:b) { 
    if(a==i){
      alpha_star[a]<-((n_t+b^2)*(b+n_new-1)*yB...[i]-(b+n_new-n_t-b^2)*y.M..[j]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
    }
    if(a!=i){
      alpha_star[a]<-(yB...[a])/(b+n_new)-mu_star
    }
  }
  for (d in 1:b) { 
    if(d==j){
      beta_star[d]<-((n_t+b^2)*(b+n_new-1)*y.M..[j]-(b+n_new-n_t-b^2)*yB...[i]-(b+n_new)*(b+n_new-1)*y....)/((b+n_new)*((b+n_new)+(n_t+b^2)*(b+n_new-2)))
    }
    if(d!=j){
      beta_star[d]<-(y.M..[d])/(b+n_new)-mu_star
    }
  } 
  R_reduce=y....*mu_star+sum(yB...*alpha_star)+sum(y.M..*beta_star)
  rowSSunadj1=sum(yB...^2/k1)-(y....)^2/N
  columnSSunadj1=sum(y.M..^2/k2)-(y....^2/N)
  totalSS1=t(y)%*%y-(y....^2/N)
  treatmentSSunadj=sum(yB..k^2/k3)+sum(yBMjg^2)-y....^2/(N)
  treatmentSSadj=R_full- R_reduce
  errorSS=totalSS1- rowSSunadj1-columnSSunadj1-treatmentSSadj
  checkSS1=sum(yB..k^2/k3)-sumcheck^2/n_c
  newSS=sum(yBMjg^2)-(sum(yBMjg))^2/n_t
  New_vs._checkSS=treatmentSSunadj-checkSS1-newSS
  dfE=N-2*b-m+2 
  meanT=treatmentSSadj/(m-1)
  meanE=errorSS/dfE
  meanCheck1=checkSS1/(c-1)
  meanNew=newSS/(n_t-1)
  meanNew_vs._checkSS=New_vs._checkSS/(1)
  treatmentF=meanT/meanE
  checkF=meanCheck1/meanE
  newF=meanNew/meanE
  New_vs._checkSSF=meanNew_vs._checkSS/meanE
  p.value.treatment=pf(treatmentF, df1=(m-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.check=pf(checkF, df1=(c-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new=pf(newF, df1=(n_t-1), df2=dfE, ncp=0, lower.tail = FALSE)
  p.value.new_vs._check=pf(New_vs._checkSSF, df1=(m-c), df2=dfE, ncp=0, lower.tail = FALSE)
  df=c(m-1,c-1,n_t-1,1,dfE)
  ss=c(treatmentSSunadj,checkSS1,newSS,New_vs._checkSS,errorSS)
  ms=c(meanT,meanCheck1,meanNew,meanNew_vs._checkSS,meanE)
  F=c(NA,NA,newF,New_vs._checkSSF,NA)
  p.value=c(NA,NA,p.value.new,p.value.new_vs._check,NA)
  anova.table=cbind(df,ss,ms,F, p.value)
  rownames(anova.table)=c("Treatment(Unadj.)","  Check(Unadj.)","  New " ,"  New vs. check","Error")
  colnames(anova.table)=c("Df", "Sum Sq","Mean Sq", "F-value","P-value")
  return(anova.table)
} 
###########
###########
###########
aov.AILSD = function(obs,row,column,treatment,type, y,n_row,n_column,n_check) 
{
  output1=aov.AILSD1(obs,row,column,treatment,type, y,n_row,n_column,n_check)  
  output2=aov.AILSD2(obs,row,column,treatment,type, y,n_row,n_column,n_check) 
  newlist <- list(output1,output2) 
  return(newlist)
}
