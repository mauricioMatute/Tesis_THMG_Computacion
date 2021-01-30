
install.packages("corrplot")
install.packages("factoextra")

x <- read.csv("base_train1_prueba_2.csv", header=TRUE)

head(x)
summary(x)
base2<-x[,c(1:9,12)]

table(base2$especie)

library(corrplot)
corrplot.mixed(cor(base2[,1:9]), order="hclust", tl.col="black")
base_bac<- base2[base2$especie == 0,]
base_hum<- base2[base2$especie == 1,]

corrplot.mixed(cor(base_bac[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_hum[,1:9]), order="hclust", tl.col="black")


dat.pca1 <- prcomp(base2[,c(1:9)], center = TRUE,scale. = TRUE)


summary(dat.pca1)
summary(base2)
sd(base2$GC_content)
sd(base2$CpG_ratio)
sd(base2$dws)
sd(base2$dry)
sd(base2$dmk)
sd(base2$LTP)
sd(base2$HTP)
sd(base2$VTP)
sd(base2$ITP)


dat.pca1$rotation
matri<-dat.pca1$x

library(ggplot2)
library(factoextra)
fviz_eig(dat.pca1)

fviz_pca_var(dat.pca1, axes=c(1,3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(dat.pca1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

