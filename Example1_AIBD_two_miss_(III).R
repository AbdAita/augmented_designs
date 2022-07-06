# Example1:An augmented randomized incomplete block design using real data (ARCBD with two missing checks safflower experiment)#
block<-c(rep(1,33),rep(2,33),rep(3,33))
treatment<-c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12","G13","G15","G16","G17","G18","G20","G21","G22","G23","G24","G25","G26","G27","G14","G19","G35","G47","G58","G83","G84","G85","G33","G34","G36","G37","G38","G39","G40","G41","G42","G43","G44","G45","G46","G48","G49","G50","G51","G52","G53","G54","G55","G56","G57","G59","G60","G14","G19","G35","G47","G58","G83","G84","G85","G67","G68","G69","G70","G71","G72","G73","G74","G75","G76","G77","G78","G79","G80","G81","G82","G86","G87","G88","G89","G90","G91","G92","G93","G94","G14","G19","G35","G47","G58","G83","G84","G85")
type<-c("test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","check","check","check","check","check","check","check","check","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","check","check","check","check","check","check","check","check","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","test","check","check","check","check","check","check","check","check")
y<-c(442.5,337.1,503.4,680,739.35,655,605,455.5,455,958.05,740.4,374.2,665.7,661.55,749.25,532.45,870.85,654.05,591.65,527.25,585,527,848,446.4,899.45,538,1001.9,856.5,384.65,971.75,800.75,469.15,670,1275.6,720,507.7,597,290.35,499.4,335.5,695,391.55,168.35,325.55,435.5,530.3,520.6,333.95,601.8,731,894.25,313.2,398,408.7,497.5,180.7,596.25,700,392.6,724.1,850,533.4,939.5,870,1037.05,805.75,434.25,775,539.5,784.5,378.15,603.1,564,596.6,919.05,490.75,688.45,751.35,643.05,792.75,305.65,508,545.85,562.9,422.9,444.4,795,837.5,639.7,710,585.45,484.05,505,478.5,523.95,459.6,537.05,586.15,205)
augmented_AIBD<-data.frame(block,treatment,type,y)
augmented_AIBD<-augmented_AIBD[-c(66),]#The check treatment G35 has been missed in the third block#
augmented_AIBD<-augmented_AIBD[-c(93),]#The check treatment G85 has been missed in the third block#
block<-augmented_AIBD[,1]
treatment<-augmented_AIBD[,2]
type<-augmented_AIBD[,3]
y<-augmented_AIBD[,4]
obs<-c(1:97)
augmented_AIBD<-data.frame(obs,block,treatment,type,y)
TMC.III.ARIBD(obs, block, treatment,type,y,2,3,8,3)