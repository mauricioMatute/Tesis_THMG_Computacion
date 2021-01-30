install.packages("corrplot")
setwd("/Users/paola.mateos/Desktop/Tesis_Mau")

getwd()

base <- read.csv("base.csv", header=TRUE)
summary(base)

sd(base$GC_content)
sd(base$CpG_ratio)
sd(base$dws)
sd(base$dry)
sd(base$dmk)
sd(base$LTP)
sd(base$HTP)
sd(base$ITP)
sd(base$VTP)
sd(base$especie)

base <- read.csv("baseTest4_pr.csv", header=TRUE)
base2<-base[,c(1:5,7,9,12)] #selección de variables
names(base2)
summary(base2)

#base_scale<-scale(base2, center = TRUE, scale = TRUE)
#summary(base_scale)

base3<-base2[base2$CpG_ratio<=20,] #remoción de outliers -> Generación de la BD
sd(base3$CpG_ratio)

base_scale<-scale(base3, center = c(0,0,0,0,0,0,0,0), scale = c(1,20,2,2,2,1,1,1)) #reducción -> Modelado
summary(base_scale)



base_scale_df<-as.data.frame(base_scale)
write.csv(base_scale_df, "baseTest4_prep.csv",row.names = FALSE)

base_CpG<-base2[base2$CpG_ratio>=20,]
summary(base_CpG)