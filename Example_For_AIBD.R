# Example1:This example has been taken from Federer(1956)#
block<-c(rep(1,7),rep(2,6),rep(3,7))
treatment<-c("H8","C3","C4","H3","C1","C2","H7","C4","C2","C1","C3","H1","H5","H4","C3","C1","H2","C4","C2","H6")
type<-c("test","check","check","test","check","check","test","check","check","check","check","test","test","test","check","check","test","check","check","test")
y<-c(74,78,78,70,83,77,75,91,81,79,81,79,78,96,87,92,89,81,79,82)
augmented_AIBD<-data.frame(block,treatment,type,y)
augmented_AIBD
augmented_AIBD<-augmented_AIBD[-c(18),]#For example,the missing value in the third block and check treatment C4#
block<-augmented_AIBD[,1]
treatment<-augmented_AIBD[,2]
type<-augmented_AIBD[,3]
y<-augmented_AIBD[,4]
obs<-c(1:19)
augmented_AIBD<-data.frame(obs,block,treatment,type,y)
augmented_AIBD
aov.ARIBD(obs, block, treatment,type,y,3,4)
