# correlation between variables
library("Hmisc")
library(corrplot)
varmatrix <- read_tsv("WorldClim_table.tsv") %>% column_to_rownames(., var="pop")
snowdata <- read_tsv("Snow_table.tsv") %>% column_to_rownames(., var="variable")
varmatrix <- rbind(varmatrix, snowdata)

varmatrix2 <- t(as.matrix(varmatrix))

varlist <- varmatrix %>% pivot_longer(cols = 2:9, names_to = "pop")
  
corplo1 <- (rcorr(as.matrix(varmatrix2)))
corplo1$P[] < 0.05 & corplo1$r[]>0.3

corplo2 <- cor(x=varmatrix2)
corplo3 <- cor.mtest(varmatrix2, conf.level = 0.95)

corrplot(corplo2, p.mat = corplo3$p, type = "upper", sig.level = 0.5)
