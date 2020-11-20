# correlation between variables
library(tidyverse)
library("Hmisc")
library(corrplot)
library(viridis)
varmatrix <- read_tsv("WorldClim_table.tsv") %>% column_to_rownames(., var="pop")
snowdata <- read_tsv("Snow_table.tsv") %>% column_to_rownames(., var="variable")
varmatrix <- rbind(varmatrix, snowdata)

varmatrix2 <- t(as.matrix(varmatrix))

corplo1 <- (rcorr(as.matrix(varmatrix2)))
corplo1$P[] < 0.05 & corplo1$r[]>0.3

corplo2 <- cor(x=varmatrix2)
corplo3 <- cor.mtest(varmatrix2, conf.level = 0.95)

corrplot(corplo2, p.mat = corplo3$p, type = "upper", sig.level = 0.5)

clust2 <- hclust(as.dist(1 - corplo2), method = "ward.D")
membs <- cutree(clust2, k=6)
plot(clust, hang = -1, cex = 0.6, ylab = NULL)


# PCA of variables
varpca <- prcomp(varmatrix2, scale. = TRUE) #SE

summary(varpca)
biplot(varpca, xlab = "PC1: 54.85% of variance", ylab = "PC2: 19.3% of variance", col = c("navy", "deeppink3"), pc.biplot = FALSE)
plotcoords=data.frame(varpca$x)
envloading1=data.frame(varpca$rotation)
envloading <- envloading1*10
PCAvar <- ggplot()+
  geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) ,y=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                   xend= c(envloading[1,1],envloading[2,1],envloading[3,1],envloading[4,1],envloading[5,1],envloading[6,1],envloading[7,1],envloading[8,1],envloading[9,1],envloading[10,1],envloading[11,1],envloading[12,1],envloading[13,1],envloading[14,1], envloading[15,1],envloading[16,1], envloading[17,1], envloading[18,1], envloading[19,1], envloading[20,1], envloading[21,1]), 
                   yend=c(envloading[1,2],envloading[2,2],envloading[3,2],envloading[4,2],envloading[5,2],envloading[6,2],envloading[7,2],envloading[8,2],envloading[9,2],envloading[10,2],envloading[11,2],envloading[12,2],envloading[13,2],envloading[14,2], envloading[15,2],envloading[16,2], envloading[17,2], envloading[18,2], envloading[19,2], envloading[20,2], envloading[21,2])), 
               color="#999999", arrow=arrow(angle = 20, length = unit(0.1,"cm"), ends = "last", type = "open"), size = 1)+
  scale_color_viridis()+
  #geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC2, color = varmatrix2[,1], label=rownames(varmatrix2)), size = 4)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 55% of variance")+
  ylab("PC 2 -- 19% of variance")+
  labs(colour = "bio1")+
  annotate("text", 
           x=c(envloading[1,1]-0.8,
               envloading[2,1],
               envloading[3,1]-0.15,
               envloading[4,1],
               envloading[5,1]-0.7,
               envloading[6,1],
               envloading[7,1]+0.25,
               envloading[8,1],
               envloading[9,1],
               envloading[10,1],
               envloading[11,1],
               envloading[12,1],
               envloading[13,1],
               envloading[14,1]-1,
               envloading[15,1],
               envloading[16,1]-0.25,
               envloading[17,1],
               envloading[18,1]+0.15,
               envloading[19,1],
               envloading[20,1]+0.5,
               envloading[21,1]), 
           y=c(envloading[1,2]+0.15,
               envloading[2,2]+0.15,
               envloading[3,2]+0.15,
               envloading[4,2]-0.15,
               envloading[5,2]-0.15,
               envloading[6,2]+0.15,
               envloading[7,2]+0.15,
               envloading[8,2]-0.15,
               envloading[9,2]+0.15,
               envloading[10,2]-0.15,
               envloading[11,2]+0.15,
               envloading[12,2]+0.15,
               envloading[13,2]+0.15,
               envloading[14,2]-0.15,
               envloading[15,2]+0.15,
               envloading[16,2]+0.15,
               envloading[17,2]-0.15,
               envloading[18,2]+0.15,
               envloading[19,2]-0.15,
               envloading[20,2]-0.15,
               envloading[21,2]+0.15), 
           label =  c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","jan_depth","snow_days"), color = "#999999")
PCAvar

# PCA of results
varpca <- prcomp(matrix1, scale. = TRUE) #SE

summary(varpca)
biplot(varpca, xlab = "PC1: 36.84% of variance", ylab = "PC2: 19.46% of variance", col = c("navy", "deeppink3"), pc.biplot = FALSE)
plotcoords=data.frame(varpca$x)
envloading1=data.frame(varpca$rotation)
envloading <- envloading1*10

PCAvar <- ggplot()+
 geom_segment(aes(x=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) ,y=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
                  xend= c(envloading[1,1],envloading[2,1],envloading[3,1],envloading[4,1],envloading[5,1],envloading[6,1],envloading[7,1],envloading[8,1],envloading[9,1],envloading[10,1],envloading[11,1],envloading[12,1],envloading[13,1],envloading[14,1], envloading[15,1],envloading[16,1], envloading[17,1], envloading[18,1], envloading[19,1], envloading[20,1], envloading[21,1]), 
                  yend=c(envloading[1,2],envloading[2,2],envloading[3,2],envloading[4,2],envloading[5,2],envloading[6,2],envloading[7,2],envloading[8,2],envloading[9,2],envloading[10,2],envloading[11,2],envloading[12,2],envloading[13,2],envloading[14,2], envloading[15,2],envloading[16,2], envloading[17,2], envloading[18,2], envloading[19,2], envloading[20,2], envloading[21,2])), 
              color="#999999", arrow=arrow(angle = 20, length = unit(0.1,"cm"), ends = "last", type = "open"), size = 1)+
  scale_color_viridis()+
  geom_point(aes(x=plotcoords$PC1, y=plotcoords$PC2, color = matrix1[,1]), size = 4)+
  theme_bw()+
  theme(axis.text=element_text(size=16))+
  theme(axis.title.x = element_text( size = 16))+
  theme(axis.title.y = element_text(size = 16))+
  xlab("PC 1 -- 36.84% of variance")+
  ylab("PC 2 -- 19.46% of variance")+
  labs(colour = "bio1")+
  annotate("text", 
           x=c(envloading[1,1],
               envloading[2,1],
               envloading[3,1],
               envloading[4,1],
               envloading[5,1],
               envloading[6,1],
               envloading[7,1],
               envloading[8,1],
               envloading[9,1],
               envloading[10,1],
               envloading[11,1],
               envloading[12,1],
               envloading[13,1],
               envloading[14,1],
               envloading[15,1],
               envloading[16,1],
               envloading[17,1],
               envloading[18,1],
               envloading[19,1],
               envloading[20,1],
               envloading[21,1]), 
           y=c(envloading[1,2],
               envloading[2,2],
               envloading[3,2],
               envloading[4,2],
               envloading[5,2],
               envloading[6,2],
               envloading[7,2],
               envloading[8,2],
               envloading[9,2],
               envloading[10,2],
               envloading[11,2],
               envloading[12,2],
               envloading[13,2],
               envloading[14,2],
               envloading[15,2],
               envloading[16,2],
               envloading[17,2],
               envloading[18,2],
               envloading[19,2],
               envloading[20,2],
               envloading[21,2]), 
           label =  c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19","jan_depth","snow_days"), color = "#999999")
PCAvar
