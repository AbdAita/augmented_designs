aov1.ARIBD = function(obs, block, treatment,type, y,W,Z) 
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
  i=W
  j=Z
  mu_b=((b*c-b-c)*sumcheck+c*(yBM..[j])+b*yBi..[i])/(b*c*(b-1)*(c-1))
  P=((sum(yBi..^2)-(yBi..[i]^2))/c)+((b*yBi..[i]-(sum(yBM..)-yBM..[j]))*(yBi..[i]+(yBM..[j])/(b-1))/(b*(c-1)))-(mu_b*(sum(yBi..)-yBi..[i]))
  TT=(sum(yBM..^2)-yBM..[j]^2)/b+yBM..[j]^2/(b-1)+sum(yBMjg^2)-sum(yB...^2/(c+h))+(yB...[i]^2/(c+h[i]))-(yB...[i]^2/(c+h[i]-1))+P
  blockSSunadj=sum(yB...^2/k1)-y....^2/(N)
  totalSS=t(y)%*%y-y....^2/(N)
  treatmentSSadj=TT 
  errorSS=totalSS- blockSSunadj-treatmentSSadj
  checkSS=sum(yBM..^2)/b+sum(yBi..^2)/c+((c*yBM..[j]+b*yBi..[i]-sumcheck)^2)/(c*b*(c-1)*(b-1))-(sumcheck^2)/(c*b)-(sum(yBi..^2)-(yBi..[i]^2))/c-(yBi..[i]^2)/(c-1)
  newSS=sum(yBMjg^2)-(sum(yBMjg))^2/n_t
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
aov2.ARIBD = function(obs, block, treatment,type, y,W,Z) 
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
  i=W
  j=Z
  mu_b=((b*c-b-c)*sumcheck+c*(yBM..[j])+b*yBi..[i])/(b*c*(b-1)*(c-1))
  P=((sum(yBi..^2)-(yBi..[i]^2))/c)+((b*yBi..[i]-(sum(yBM..)-yBM..[j]))*(yBi..[i]+(yBM..[j])/(b-1))/(b*(c-1)))-(mu_b*(sum(yBi..)-yBi..[i]))
  TT=(sum(yBM..^2)-yBM..[j]^2)/b+yBM..[j]^2/(b-1)+sum(yBMjg^2)-sum(yB...^2/(c+h))+(yB...[i]^2/(c+h[i]))-(yB...[i]^2/(c+h[i]-1))+P
  treatmentSSunadj=sum(y.M..^2/k)-y....^2/(N)
  totalSS=t(y)%*%y-y....^2/(N)
  blockSSadj=P
  errorSS=totalSS-blockSSadj-treatmentSSunadj
  checkSS1=(sum(yBM..^2)-yBM..[j]^2)/b+yBM..[j]^2/(b-1)-(sumcheck^2/n_c)
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
aov.ARIBD = function(obs, block, treatment,type, y, W, Z)
{
  output1=aov1.ARIBD(obs, block, treatment,type, y, W, Z) 
  output2=aov2.ARIBD(obs, block, treatment,type, y, W, Z)
  newlist <- list(output1,output2) 
  return(newlist)
}

