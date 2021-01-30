
install.packages("corrplot")
install.packages("factoextra")

base <- read.csv("base.csv", header=TRUE)

head(base)
base2<-base[,c(1:9,12)]

library(corrplot)
corrplot.mixed(cor(base2[,1:9]), order="hclust", tl.col="black")
base_lev <- base2[base2$especie == 0,]
base_ciona <- base2[base2$especie == 1,]
base_mosca <- base2[base2$especie == 2,]
base_pezC <- base2[base2$especie == 3,]
base_pollo <- base2[base2$especie == 4,]
base_rat <- base2[base2$especie == 5,]
base_hum <- base2[base2$especie == 6,]
base_gus <- base2[base2$especie == 7,]

corrplot.mixed(cor(base_lev[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_ciona[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_mosca[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_pezC[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_pollo[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_rat[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_hum[,1:9]), order="hclust", tl.col="black")
corrplot.mixed(cor(base_gus[,1:9]), order="hclust", tl.col="black")


dat.pca1 <- prcomp(base2[,c(1:9)], center = TRUE,scale. = TRUE)
x <-summary(base2)

sumall <- sapply(base2, function(x) c(summary(x), type = class(x)))
write.csv(sumall, file = 'sumall.csv')




summary(dat.pca1)
dat.pca1$rotation
matri<-dat.pca1$x

                 
library(ggplot2)
library(factoextra)
fviz_eig(dat.pca1)

fviz_pca_var(dat.pca1, #axes=c(2,3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


fviz_pca_var(dat.pca1, axes=c(1,3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(dat.pca1, axes=c(2,3),
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
#fviz_pca_ind(dat.pca1,
#             col.ind = "cos2", # Color by the quality of representation
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE     # Avoid text overlapping
#)


