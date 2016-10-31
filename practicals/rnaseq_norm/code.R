
# First we need to install Bioconductor
source("http://bioconductor.org/biocLite.R")


# We need the DESeq library
if(!required(DESeq)){
	library(DESeq)
}

# We need the marray library
if(!required(marray)){
  biocLite()
}


# ## Location of datasets
# url.pheno <- "ftp://tagc.univ-mrs.fr/upload/parathyroid_pheno.txt"
# url.counts <- "ftp://tagc.univ-mrs.fr/upload/parathyroidGenes.txt"
# 
# ## Loading count table
# cnt <-  read.table(url.counts,sep="\t", head=T, row=1)
# #dimensions
# dim(cnt)
# #first 2 lines
# head(cnt,2)
# 
# ## Load phenotypic data
# pheno <- read.table(url.pheno, sep='\t', head=TRUE, row=1)
# 
# 
# 
# 
# dir.output <- "~/Parathyroid" ## Define the output directory. You can adapt this to your local configuration.
# dir.create(dir.output, showWarnings=FALSE, recurs=TRUE) ## Create the directory if it does not exist yet
# setwd(dir.output)
# 
# 
# 
# ## Pheno object dimension
# dim(pheno)
# ncol(pheno)
# nrow(pheno)
# 
# ## Yes all samples are referenced.
# all(colnames(cnt) == rownames(pheno))
# 
# ## Changing colum names
# colnames(cnt) <- paste(
#   colnames(cnt), 
#   pheno$patient, 
#   pheno$time, 
#   pheno$treatment, sep="|"
# )
# colnames(cnt)
# 
# # Deleting genes with 0 counts
# row.sum <- apply(cnt,1,sum)
# cnt <- cnt[ row.sum > 0, ]
# 
# # Extracting Control/DPN treament at 24h.
# 
# # 
# # control <- cnt[, pheno$time=="24h" & pheno$treatment=="Control"]
# # dpn     <- cnt[, pheno$time=="24h" & pheno$treatment=="DPN"]
# # cnt.sub <- cbind(control, dpn)
# # dim(cnt.sub)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## Now we can start by loading the data.set of interest 
# ## The read.table function can be used alternatively 
# ## if the dataset comes as a flat file.
# 
# load(DESeq)
# load(parathyroid)
# 
# ## Let's get some informations about the parathyroid package:
# ## NB: q to quit the description.
# help(package="parathyroid")
#  
# ## We will use the parathyroidGenes dataset
# 
# data(parathyroidGenes)
# # the parathyroidGenes is an object/instance of class CountDataSet.
# is(parathyroidGenes)
# 
# slotNames(parathyroidGenes)
# help("CountDataSet-class")
# 
# # access to the counts is provided through the 'counts' method:
# head(counts(parathyroidGenes))
# counts <- counts(parathyroidGenes)
# 
# # access to the phenotypes informations is provided though phenoData
# # which returns a kind of data.frame (in fact an object of class 'AnnotatedDataFrame')
# phenoData(parathyroidGenes)
# phenotypes <- pData(phenoData(parathyroidGenes))
# colnames(phenotypes)
# #[1] "sizeFactor" "experiment" "patient"    "treatment"  "time"      
# #[6] "submission" "study"      "sample"     "run"  
# phenotypes$treatment
# #[1] Control Control DPN     DPN     OHT     OHT     Control Control DPN    
# #[10] DPN     DPN     OHT     OHT     OHT     Control Control DPN     DPN    
# #[19] OHT     OHT     Control DPN     DPN     DPN     OHT     OHT     OHT    
# #Levels: Control DPN OHT
# phenotypes$time
# #[1] 24h 48h 24h 48h 24h 48h 24h 48h 24h 24h 48h 24h 24h 48h 24h 48h 24h 48h 24h
# #[20] 48h 48h 24h 48h 48h 24h 48h 48h
# #Levels: 24h 48h
# 
# 
# 
# # Let's extract a subset of data : control vs DPN at time point 24h
# 
# control <- counts[, phenotypes$time=="24h" & treatment=="Control"]
# dpn <- counts[, phenotypes$time=="24h" & treatment=="DPN"]
# counts.controlVsDPN.24h <- cbind(control,dpn)
# dim(counts.controlVsDPN.24h)
# pheno.controlVsDPN.24h <- c(rep(0,4),rep(1,4))
# 
# # let's have a look at the data
# # A pseudo count is added to avoid warnings 
# # because of 0 values.
# 
# counts.controlVsDPN.24h.log <- log2(counts.controlVsDPN.24h+1)
# 
# #scatter plot
# plot(	 counts.controlVsDPN.24h.log[,1], 
# 	 counts.controlVsDPN.24h.log[,5], 
# 	 pch=".", 
# 	 xlab="Control_1",  
# 	 ylab="DPN_1"
# )
# 
# # distribution
# library(affy)
# 
# plotDensity(counts.controlVsDPN.24h.log,
# 	col=pheno.controlVsDPN.24h+1,
# 	lty=pheno.controlVsDPN.24h+1
# )
# 
# 
# 
# parathyroidGenes <- estimateSizeFactors( parathyroidGenes )
# is(parathyroidGenes)
# counts.n <- counts( parathyroidGenes , normalized=TRUE )
# boxplot(log2(counts.n+1), pch=".")


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


## We need the NOISeq library
# contains the data set we will play with.
if(!require(NOISeq)){
  biocLite("NOISeq")
}

## Load the libraries
library(edgeR)
library(marray)
library(NOISeq)

## Load the dataset
data(Marioni)
head(mycounts) 

## mycounts will only contain genes that have a count in 
# at least one condition.
allzero <- rowSums(mycounts>0) == 0
mycounts <- mycounts[!allzero,]


## lib.size will contain library size:
lib.size <- colSums(mycounts)


## Exploratory data analysis
# Histogram of the whole dataset.
# Note that a pseudo-count is added to each cell of the matrix.
# Data are transformed in logarithmic-scale for visualization
plot(density(as.matrix(log2(mycounts+1))))

# Histogram of corresponding to each sample
 x.density <- apply(, 2, density)
 all.x <- do.call(cbind, lapply(x.density, function(x) x$x))
    all.y <- do.call(cbind, lapply(x.density, function(x) x$y))
    matplot(all.x, all.y, ylab = ylab, xlab = xlab, type = type, 
        col = col, ...)


## upperQuartile normalization as computed in edgeR
# Note that the implementation also use the lib.size
# 1. Normalize by library size
y <- sweep(mycounts, 2, lib.size, FUN="/")
# 2. Compute upper quartile (q75)
f <- apply(y,2,quantile, 0.75)
# 3. Scaling factors should multiply to one
f.q3 <- f/exp(mean(log(f)))
# 4. Check the result is conform
f.q3 == calcNormFactors(mycounts, method="upperquartile")



# 

## 3rd quartile (Q3)

mycounts.col.q3 <- apply(mycounts,2, quantile, 0.75)
scaleFactor <- mycounts.col.q3/min(mycounts.col.q3)
mycounts.n.q3    <- sweep(  mycounts.raw, 
                              MARGIN=2, 
                              scaleFactor, 
                              FUN="/"
                      )


## Quantile normalization (Q)
mycounts.n.quant <- normalizeQuantiles(mycounts)

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

