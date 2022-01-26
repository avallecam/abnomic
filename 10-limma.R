library(tidyverse)

# kasper tutorial ---------------------------------------------------------
# https://kasperdanielhansen.github.io/genbioconductor/html/limma.html

library(leukemiasEset)
data(leukemiasEset)
leukemiasEset
biobroom::tidy.ExpressionSet(leukemiasEset,addPheno = T)

ourData <- leukemiasEset[, leukemiasEset$LeukemiaType %in% c("ALL", "NoL")]
ourData$LeukemiaType <- factor(ourData$LeukemiaType)

biobroom::tidy.ExpressionSet(ourData,addPheno = T)

biobroom::tidy.ExpressionSet(ourData,addPheno = T) %>% 
  count(LeukemiaType)

biobroom::tidy.ExpressionSet(ourData,addPheno = T) %>% 
  filter(gene=="ENSG00000163751") %>% 
  ggplot(aes(x = value, fill = LeukemiaType)) +
  geom_histogram()

library(limma)

# one exposure ------------------------------------------------------------

design <- model.matrix(~ LeukemiaType,data = ourData)
colnames(design)
fit <- 
  design %>% 
  lmFit(object = ourData, design = .) %>% 
  eBayes()
topTable(fit)

# biobroom::tidy.MArrayLM(fit)
biobroom::tidy.MArrayLM(fit) %>% arrange(p.value)
biobroom::tidy.MArrayLM(fit) %>% count(term)


# exposure and confounding ------------------------------------------------
# https://www.biostars.org/p/231771/

biobroom::tidy.ExpressionSet(leukemiasEset, addPheno = T) %>% 
  count(LeukemiaType,Subtype)
# biobroom::tidy.ExpressionSet(ourData,addPheno = T) %>%
#   count(LeukemiaType,Subtype)

design2 <- model.matrix(~ LeukemiaType + Subtype,data = leukemiasEset)
colnames(design2)
fit2 <- 
  design2 %>% 
  lmFit(object = leukemiasEset, design = .) %>%
  eBayes()

biobroom::tidy.MArrayLM(fit2)
biobroom::tidy.MArrayLM(fit2) %>% 
  count(term)

topTable(fit2,coef = "LeukemiaTypeNoL")
topTable(fit2,coef = "SubtypeT_ALL")

# smyth tutorial ----------------------------------------------------------
# https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

# URL <- "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-33005"
# SDRF.file <- "E-GEOD-33005.sdrf.txt"
# Data.file <- "E-GEOD-33005.raw.1.zip"
# download.file(paste(URL,SDRF.file,sep="/"), SDRF.file)
# download.file(paste(URL,Data.file,sep="/"), Data.file)
# unzip(Data.file)

# targets <- limma::readTargets("targets.txt")

# Batch <- factor(c(1,3,4,1,3,4))
# Pasilla <- factor(GEO$Pasilla,levels=c("Normal","Down"))
# design <- model.matrix(~ Batch + Pasilla)



# higher criticism threshold ----------------------------------------------
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6164648/

# defines a bound of random variation

# benjamini hochberg
qbinom(p = 0.95, size = 201,prob = 0.05)
qbinom(p = 0.95, size = 394,prob = 0.05)
qbinom(p = 0.95, size = 23332,prob = 0.05)

# bonferroni
qbinom(p = 1-(0.05/20), size = 23332,prob = 0.05)
qbinom(p = 1-(0.05/20), size = 394,prob = 0.05)

#kolmogorov smirnoff
x<-runif(1000,0,1)
ks.test(x,"punif",0,1)
enframe(x) %>% 
  ggplot(aes(x = value)) +
  geom_histogram()


# storey qvalue -----------------------------------------------------------


# _good histogram ---------------------------------------------------------

library(qvalue)

data(hedenfalk)
hedenfalk %>% 
  as_tibble()

pvalues <- hedenfalk$p

qobj <- qvalue(p = pvalues)
qobj %>% str()

# quick_sumviz, eval=FALSE-
summary(qobj)
hist(qobj)
plot(qobj)


# _bad histogram ----------------------------------------------------------

# pvalue_histBad, dependson=c("load_qvalue", "quick_p"), echo=FALSE, fig.height=3, fig.width=5-
set.seed(478)
p2 = c(hedenfalk$p, (runif(450, min=0.7, max=1))^(0.33))
somethingsWrong = list(p=p2)
hist(somethingsWrong$p, nclass=20, main="Problematic p-values", xlab="intentionally bad, simulated p-values")

wrong_pvalue <- somethingsWrong$p
qobj <- qvalue(p = wrong_pvalue)
summary(qobj)
hist(qobj)
plot(qobj)

