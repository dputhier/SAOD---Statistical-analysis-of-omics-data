################################################################
## handling RNA-seq data: exploratory data analysis and normalization 
##
## Original data published in
##     Marioni et al. RNA-seq: an assessment of technical
##     reproducibility and comparison with gene expression
##     arrays. Genome Res (2008) vol. 18 (9) pp. 1509-17
##

## First we need to install Bioconductor
source("http://bioconductor.org/biocLite.R")

if(!require(Biobase)){
  biocLite()
}


## We need the edgeR library
# contains some normalization methods
if(!require(edgeR)){
  biocLite("edgeR")
}

## We need the marray library
# contains the quantile normalization method
if(!require(marray)){
  biocLite("marray")
}


## Load the libraries
library(edgeR)
library(marray)

####################################################
## Load the dataset
####################################################

## Load the supplementary table 2 from Marioni et al (2008), with the
## raw counts.

url <- "http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics/data/NGS/Marioni_2008/SupplementaryTable2.tab"

marioni <- read.table(url, sep='\t', header=TRUE, row.names=1)
dim(marioni)
head(marioni)
counts <- marioni[,6:19]
 

####################################################
## pre-processing
####################################################

## counts will only contain genes that have a count in 
# at least one condition.
allzero <- rowSums(counts>0) == 0
counts <- counts[!allzero,]
dim(counts)

## What is the number of counts in each library (library size)
lib.size <- colSums(counts)
lib.size

####################################################
## Exploratory data analysis
####################################################

## Summary of each sample
# Note the maximum values that reflects the presence of high-count genes.
summary(counts)


## Distributions.
# It is generally advised two assess the distribution using with log-transformed data. 
# Log-transformed values gives us a better perception of
# the distribution, because it emphasizes the lower values. Here we will add a pseudo-count to every cell of the
# count table to avoid any warning during log-transformation

counts.log <- log2(counts+1)

## Distribution of the whole dataset
hist(as.matrix(counts.log), border="white", col="blue")
grid(col="black")
plot(density(as.matrix(counts.log)))

## Distributions for individual samples.

 mycounts.density <- apply( counts.log, 
                            MARGIN=2, 
                            FUN=density
                      )
 x.values <- do.call(cbind, lapply( mycounts.density, 
                                    function(x) x$x
                            )
                    )
y.values <- do.call(cbind, lapply( mycounts.density, 
                                    function(x) x$y
                            )
                    )
# Searching for sample class
samples.kidney <- grep("Kidney",colnames(counts))
samples.liver  <- grep("Liver",colnames(counts))
color.palette <- rep(2,ncol(counts))
color.palette[samples.kidney] <- 3
# 
matplot(x.values, y.values, type="l", col=color.palette)
grid(col="black")

##Q°:
# What can we say about distribution of counts in kidney and liver sample. What about low-count genes and high-count genes ?

# Not the same number of highly expressed genes in both conditions. -> See also reproducibility...

## Cumulative sum of counts for sample 1
counts.sort.cum <- cumsum(sort(counts[,1]))
plot(counts.sort.cum, type="l", ylab="Cumulative sum of counts", xlab="Genes ordered by count values")
grid(col="black")

# What can we guess from this graphic.

# How should we interpret this graphic. How many counts are associated with the 5% most expressed genes. 
length(which(counts.sort.cum )/lib.size[1] >= 0.5))/nrow(counts)

# less than 5% of the genes account for approximately 50% of the total counts !
# These genes may vary across biological conditions. Given a constant number of reads/counts per library, an increase in their abundance will be accompanied by a strong impact on weakly expressed genes representation.
length(which(counts.sort.cum/lib.size[1] >= 0.5))/nrow(counts)


## Comparison of one liver and one kidney sample.
# One can compare one sample from both condition using a simple scatter-plot.

plot(counts[,samples.kidney[1]], counts[,samples.kidney[2]], pch=".", col="seagreen")
grid(lwd=1, col="#000000")

# Again the use of log-transformed is well-suited here.
plot(counts.log[,samples.kidney[1]], 
      counts.log[,samples.liver[2]], 
      pch=".",
      cex=2, 
      col="seagreen", 
      xlab="Gene counts (Kidney)", 
      ylab="Gene counts (liver)"
)
grid(lwd=1, col="#000000")

## Same comparison but the density of points is color coded.
# The 0 values appear clearly in the bottom left corner.
library(geneplotter)
smoothScatter(counts.log[,samples.kidney[1]], 
              counts.log[,samples.liver[2]], 
              pch=".",
              xlab="Gene counts (Kidney)", 
              ylab="Gene counts (liver)",
              cex=2)

## It is most generally  good idea to perform a more general comparison.
# Here scatter-plot for all pair of kidney replicates are displayed:

myDisplayFunction <- function(x,y){smoothScatter(x,y,pch=".", add=TRUE)}
pairs(counts.log[,samples.kidney], 
      upper.panel=myDisplayFunction,  
      lower.panel=NULL
)

## The MA-plot is another popular way of displaying the count-results of two samples. 
# When more than two samples are available, one can use a pseudo sample as a reference. This pseudo sample here 
# contains representative values computed as the median of each gene. 

# MA plot for the first sample
pseudo.sample  <- apply(counts.log, 1, median) 
M <- counts.log - pseudo.sample
A <- (counts.log + pseudo.sample)/2

smoothScatter(A[,1],M[,1], pch=".")
abline(h=0)
abline(h=-1)
abline(h=1)
lines(lowess(M[,1]~A[,1]), col="red")

# MA plot for all the samples
par(mfrow=c(4,4))
for(i in 1:ncol(counts.log)){
     smoothScatter(A[,i],M[,i], pch=".")
     abline(h=0)
     abline(h=-1)
     abline(h=1)
     lines(lowess(M[,i]~A[,i]), col="red")
}

##M-values distributions. 
# We can check the M-value distribution that is expected to be centered on zero as expression levels of most of the genes are expected to be unchanged.
# We can also check the expression levels of some House-keeping genes.
# First we load a list of housekeeping genes.

hkg <- read.table("/home/denis/DOCUMENTS/GIT/statistics_bioinformatics/practicals/ASG1_2012/rnaseq_norm/pone.0000898.s001.C1.C2.500top.tsv", sep="\t", head=T)
hkg.EnsId <-  as.character(hkg[,2])
names(hkg.EnsId) <- hkg[,1]


# then we draw the histogram of M values. Here for the comparison of the first liver and the first kidney samples.

isHouseKeeping <- rownames(counts) %in% hkg.EnsId
br <- seq(min(M), max(M), length.out=100)
myGreen <- rgb(10,100,10,max=255, alpha=150)
myGray <- rgb(100,100,100,max=255, alpha=100)

hist(M[!isHouseKeeping,samples.kidney[1]], 
      col=myGray, 
      breaks=br,
      freq=FALSE
) 


hist( M[isHouseKeeping, samples.kidney[1]], 
      col=myGreen, 
      breaks=br,
      freq=FALSE,
      add=TRUE
) 

abline(v=0, col="black",lwd=2)
grid(col="black")



########################################################
## Normalization
########################################################
## Below are some examples of popular normalization methods used for RNA-Seq.

## Total count normalization:
#Total count (TC): For each sample, gene counts are divided by the corresponding total number of  mapped reads (library size)  then multiplied by the mean total count across all the samples of the dataset.



## upperQuartile normalization as computed in edgeR
# Note that the implementation also use the lib.size
# 1. Normalize by library size
y <- sweep(counts, 2, lib.size, FUN="/")
# 2. Compute upper quartile (q75)
f <- apply(y,2,quantile, 0.75)
# 3.  For symmetry, normalization factors are adjusted to multiply to 1.
#The effective library size is then the original library size
#multiplied by the scaling factor.
f.q3 <- f/exp(mean(log(f)))
prod(f.q3)
# 4. Check the result is conform
f.q3 == calcNormFactors(counts, method="upperquartile")

# 5. Normalized data
counts.n.q3 <- sweep(counts, 2, f.q3, FUN="/")


## median normalization as computed in edgeR
# Note that the implementation also use the lib.size
# 1. Normalize by library size
y <- sweep(mycounts, 2, lib.size, FUN="/")
# 2. Compute median
f <- apply(y,2,median)
# 3.  For symmetry, normalization factors are adjusted to multiply to 1.
f.q2 <- f/exp(mean(log(f)))
prod(f.q2)
# 4. Check the result is conform
f.q2 == calcNormFactors(mycounts, method="upperquartile"  )


## Quantile normalization (Q)
mycounts.n.quant <- normalizeQuantiles(mycounts)

## RLE normalization (as computed in edgeR)
gm <- exp(rowMeans(log(counts)))
f <- apply(data, 2, function(u) median((u/gm)[gm > 0]))
f <- f/lib.size

## Library size (LS)
lib.size    <- colSums(mycounts)
scaleFactor <- lib.size/min(lib.size)




mycounts.n.ls    <- sweep(  mycounts, 
                              MARGIN=2, 
                              scaleFactor, 
                              FUN="/"
                      )


## Mediane (MED):

mycounts.col.med <- apply(mycounts,2, median)
scaleFactor <- mycounts.col.med/min(mycounts.col.med)
mycounts.n.med    <- sweep(  mycounts.raw, 
                              MARGIN=2, 
                              scaleFactor, 
                              FUN="/"
                      )

## Median (computed with values > 0,  MEDnz):

mycounts.col.med <- apply(mycounts.raw.zero.to.NA,2, median, na.rm=TRUE)
scaleFactor <- mycounts.col.med/min(mycounts.col.med)
mycounts.n.med.no.0    <- sweep(  mycounts.raw, 
                              MARGIN=2, 
                              scaleFactor, 
                              FUN="/"
                      )



## 3rd quartile (computed with values > 0, Q3nz )

mycounts.col.q3 <- apply(mycounts.raw.zero.to.NA,2, quantile, 0.75, na.rm=T)
scaleFactor <- mycounts.col.q3/min(mycounts.col.q3)
mycounts.n.q3.no.0    <- sweep(  mycounts.raw, 
                              MARGIN=2, 
                              scaleFactor, 
                              FUN="/"
                      )


## TMM
library(edgeR)
scaleFactor <- calcNormFactors(mycounts.raw, method="TMM")
mycounts.n.TMM <- sweep(  mycounts.raw, 
                           MARGIN=2, 
                           scaleFactor, 
                           FUN="/"
                  )


## Let's load a set of 1000 housekeeping genes

hkg <- read.table("pone.0000898.s001.C1.C2.1000top.tsv", sep="\t", head=T)
hkg.EnsId <-  hkg[,2]
hkg.EnsId <-  as.character(hkg[,2])
names(hkg.EnsId) <- hkg[,1]


## A function that will perform some diagnosis according to 
## Normalization results

drawDiagnosisPlot <- function(x=NULL, title=NULL, pair=1){
  x <- log2(x+1)
  kidney.pos  <- grep("Kidney", colnames(x))
  liver.pos   <- grep("Liver", colnames(x))
  col.boxplot <- rep(2,ncol(x))
  col.boxplot[kidney.pos] <- 3 
  X11()
  par(mfrow=c(3,2))
  boxplot(x, 
          pch=".", 
          col=col.boxplot, 
          main=title,
          las=2,
          cex.axis=0.7
  )

  br <- seq(min(x), max(x), length.out=100)
  hist(x[,kidney.pos[pair]], col=3, breaks=br)
  hist(x[,liver.pos[pair]], col=2, breaks=br)

  M <- x[,kidney.pos[pair]] - x[,liver.pos[pair]]
  names(M) <- rownames(x)

  A <- (x[,kidney.pos[pair]] + x[,liver.pos[pair]])/2
  names(A) <- rownames(x)

  myGray <- rgb(100,100,100,max=255, alpha=100)
  myGreen <- rgb(10,100,10,max=255, alpha=255)

  plot( A,
        M,
        pch=".",
        col=myGray,
        main=paste(title, " (pair: ", pair, ")")
  )
  
  abline(h=0)

  isHouseKeeping <- rownames(x) %in% hkg.EnsId
  
  

  points( A[isHouseKeeping],
          M[isHouseKeeping],
          col=myGreen,
          pch=16
    )

  lines(  lowess(M~A), 
          col="red", 
          lwd=2
  )

  br <- seq(min(M), max(M), length.out=100)

  myGreen <- rgb(10,100,10,max=255, alpha=150)

  hist(M[!isHouseKeeping], 
        col=myGray, 
        breaks=br,
        freq=FALSE,
        main=paste(title, " (pair: ", pair, ")")
  ) 


  hist( M[isHouseKeeping], 
        col=myGreen, 
        breaks=br,
        freq=FALSE,
        add=TRUE
  ) 

  abline(v=0, col="black",lwd=2)
  grid(col="black")

}

drawDiagnosisPlot(mycounts.raw, title="Raw (log-transformed)")
drawDiagnosisPlot(mycounts.n.quant, title="Quantile normalization")
drawDiagnosisPlot(mycounts.n.med, title="Median" )
drawDiagnosisPlot(mycounts.n.med.no.0, title="Median (compute on values > 0)" )
drawDiagnosisPlot(mycounts.n.q3, title="3rd Quartile" )
drawDiagnosisPlot(mycounts.n.q3.no.0, title="3rd Quartile (compute on values > 0)" )
drawDiagnosisPlot(mycounts.n.ls, title="Library size" )
drawDiagnosisPlot(mycounts.n.TMM, title="TMM" )
