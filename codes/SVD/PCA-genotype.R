# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("SNPRelate")

library("SNPRelate")

genofile <- snpgdsOpen(snpgdsExampleFileName())

# run PCA
RV <- snpgdsPCA(genofile)
RV

# eigenvalues
head(RV$eigenval)

# variance proportion (%)
head(round(RV$varprop*100, 2))
# [1] 12.23  5.84  1.01  0.95  0.84  0.74

# draw
plot(RV)
plot(RV, 1:4)


####  there is no population information  ####

# make a data.frame
tab <- data.frame(sample.id = RV$sample.id,
                  EV1 = RV$eigenvect[,1],    # the first eigenvector
                  EV2 = RV$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")


####  there are population information  ####

# get population information
#   or pop_code <- scan("pop.txt", what=character())
#   if it is stored in a text file "pop.txt"
pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot/pop.group"))

# get sample id
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))



# assume the order of sample IDs is as the same as population codes
cbind(samp.id, pop_code)
#        samp.id       pop_code
#   [1,] "NA19152"     "YRI"   
#   [2,] "NA19139"     "YRI"   
#   [3,] "NA18912"     "YRI"   
#   [4,] "NA19160"     "YRI"   
#   [5,] "NA07034"     "CEU"   
#   ...

# make a data.frame
tab <- data.frame(sample.id = RV$sample.id,
                  pop = factor(pop_code)[match(RV$sample.id, samp.id)],
                  EV1 = RV$eigenvect[,1],    # the first eigenvector
                  EV2 = RV$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)


# draw
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop),
     xlab="eigenvector 2", ylab="eigenvector 1")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:4)


#SVD
A <- SNPRelate::snpgdsGetGeno(genofile)
A <- scale(A,center = T, scale = F)
svd(A)
A <- A[,-unique(which(is.na(A),arr.ind = T)[,2])]
svdA <- svd(A,nu = 4, nv = 4)
pairs(svdA$u)
plot(RV, 1:4)

# close the file
snpgdsClose(genofile)