install.packages("corrplot")

base <- read.csv("base_train1_prueba_2.csv", header=TRUE)
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

base <- read.csv("base_test1_prueba_2.csv", header=TRUE)
#base2<-base[,c(1:5,7,9,12)] #selección de variables train
base2<-base[,c(1:5,7,9)] #selección de variables test
names(base2)
summary(base2)

base3<-base2[base2$CpG_ratio<=3,] #remoción de outliers -> Generación de la BD
sd(base3$CpG_ratio)

#base_scale<-scale(base3, center = c(0,0,0,0,0,0,0,0), scale = c(1,3,2,2,2,1,1,1)) #reducción train-> Modelado
base_scale<-scale(base3, center = c(0,0,0,0,0,0,0), scale = c(1,3,2,2,2,1,1)) #reducción test-> Modelado
summary(base_scale)



base_scale_df<-as.data.frame(base_scale)
write.csv(base_scale_df, "base_test1_prueba_2_prep.csv",row.names = FALSE)

base_CpG<-base2[base2$CpG_ratio>=3,]
summary(base_CpG)
