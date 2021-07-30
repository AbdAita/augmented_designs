# Example2:package plantbreeding in R language (2012)#
row<-c(rep(1,15),rep(2,15),rep(3,15),rep(4,15),rep(5,15))
column<-c(rep(1:5,15))
treatment<-c("C","44","25","E","10","19","8","B","2","D","37","A","35","26","24","46","39","A","27","5","12","C","16","D","45","E","4","11","28","B","34","15","C","7","20","D","1","42","B","38","21","E","36","32","A","17","D","E","22","C","B","13","30","A","49","6","31","9","41","33","A","18","47","3","40","23","48","D","50","E","43","B","29","C","14")
type<-c("check","test","test","check","test","test","test","check","test","check","test","check","test","test","test","test","test","check","test","test","test","check","test","check","test","check","test","test","test","check","test","test","check","test","test","check","test","test","check","test","test","check","test","test","check","test","check","check","test","check","check","test","test","check","test","test","test","test","test","test","check","test","test","test","test","test","test","check","test","check","test","check","test","check","test")
y<-c(2.1,7.36,3.67,4.3,5.45,2.69,4.03,3.89,3.8,4.45,4.56,2.88,6.61,3.2,4.2,4.26,4.1,4.42,8.08,6.7,3.39,2.95,3.68,9.95,4.65,1.02,3.18,5.22,3.23,5.9,2.17,3.24,2.57,5.3,6,9.06,3.29,1.87,5.52,7.15,2.49,6.69,3.4,6.94,4.6,4.28,8.53,7.84,5.13,3.5,5.38,2.86,3.9,5.82,7.2,3.6,2.96,3.98,4.77,4.05,2.7,2.84,6.68,9.56,6.55,4.66,5.47,7.84,5.57,4.5,6.65,5.16,2.57,3.2,6.45)
augmented_AILSD<-data.frame(row,column,treatment,type,y)
augmented_AILSD<-augmented_AILSD[-c(39),] #For example,the missing value in the third row and fourth column and check treatment B#
row<-augmented_AILSD[,1]
column<-augmented_AILSD[,2]
treatment<-augmented_AILSD[,3]
type<-augmented_AILSD[,4]
y<-augmented_AILSD[,5]
obs<-c(1:74)
augmented_AILSD<-data.frame(obs,row,column,treatment,type,y)
augmented_AILSD
aov.AILSD(obs,row,column,treatment,type, y,n_row=3,n_column=4,n_check=2 )  
