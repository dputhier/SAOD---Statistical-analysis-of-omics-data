calcNormFactors <- function(object, method=c("TMM","RLE","upperquartile","none"), refColumn=NULL, logratioTrim=.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
#       Scale normalization of RNA-Seq data.
#       Mark Robinson.  Edits by Gordon Smyth.
#       Created October 22 October 2009.  Last modified 16 Apr 2013.
{
#       Check object
        if(is(object,"DGEList")) {
                x <- as.matrix(object$counts)
                lib.size <- object$samples$lib.size
        } else {
                x <- as.matrix(object)
                lib.size <- colSums(x)
        }

#       Check method    
        method <- match.arg(method)

#       Remove all zero rows
        allzero <- rowSums(x>0) == 0
        if(any(allzero)) x <- x[!allzero,,drop=FALSE]

#       Degenerate cases
        if(nrow(x)==0 || ncol(x)==1) method="none"

#       Calculate factors
        f <- switch(method,
                TMM = {
                        f75 <- calcFactorQuantile(data=x, lib.size=lib.size, p=0.75) ## calculate norm fact with upperQ method.
                        if( is.null(refColumn) )
                          refColumn <- which.min(abs(f75-mean(f75))) ## find the sample which is the most representative (whose q75 is the closest to the mean value of q75)
                        if(length(refColumn)==0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
                        f <- rep(NA,ncol(x))   
                        for(i in 1:ncol(x))  ## for each sample
                          f[i] <- .calcFactorWeighted(	obs=x[,i],
										ref=x[,refColumn], 
										libsize.obs=lib.size[i], 
										libsize.ref=lib.size[refColumn], ## the scaling factor will be computed against this reference.
										logratioTrim=logratioTrim, 
										sumTrim=sumTrim, 
										doWeighting=doWeighting, 
										Acutoff=Acutoff)
                        f
                },
                RLE = .calcFactorRLE(x)/lib.size,
                upperquartile = .calcFactorQuantile(x,lib.size,p=p),
                none = rep(1,ncol(x))
        )

#       Factors should multiple to one
        f <- f/exp(mean(log(f)))

#       Output
        if(is(object, "DGEList")) {
                object$samples$norm.factors <- f
                return(object)
        } else {
                return(f)
        }
}


.calcFactorRLE <- function (data)
{
        gm <- exp(rowMeans(log(data)))
        apply(data, 2, function(u) median((u/gm)[gm > 0]))
}

.calcFactorQuantile <- function (data, lib.size, p=0.75)
{
#       i <- apply(data<=0,1,all)
#       if(any(i)) data <- data[!i,,drop=FALSE]
        y <- t(t(data)/lib.size)
        f <- apply(y,2,function(x) quantile(x,p=p))
#       f/exp(mean(log(f)))
}

.calcFactorWeighted <- function(	obs, ## the counts in the current sample 
							ref,  ## the counts in the reference sample
							libsize.obs=NULL, ## lib size 
							libsize.ref=NULL,  ## lib size of the reference
							logratioTrim=.3,    ## amount of trim to use on log-ratios ("M" values). We will try to focus on non differentially expressed genes
							sumTrim=0.05, ## amount of trim to use on the combined absolute levels ("A" values). Genes too low or high won't be taken into account. 
							doWeighting=TRUE, 
							Acutoff=-1e10)
{
        if( all(obs==ref) ) return(1) ## if the current sample is the reference, don't compute the scaling factor is 1.

        obs <- as.numeric(obs) ## data should already be numeric...
        ref <- as.numeric(ref)  ## data should already be numeric...

        if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs ## n0 is the libsize of current sample
        if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref     ## nR is the libsize of current sample

        logR <- log2((obs/nO)/(ref/nR))                 ## log ratio of expression (M), accounting for library size
        absE <- (log2(obs/nO) + log2(ref/nR))/2  ## absolute expression   (A values)
        v <- (nO-obs)/nO/obs + (nR-ref)/nR/ref      ## estimated asymptotic variance
									## 	
#       remove infinite values, cutoff based on A
        fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff) ## Log 0 gives infinite values (should be excluded)

        logR <- logR[fin]  ## exclude infinite values
        absE <- absE[fin] ## exclude infinite values
        v <- v[fin]  ## exclude infinite values

#       taken from the original mean() function

        n <- sum(fin)
        loL <- floor(n * logratioTrim) + 1
        hiL <- n + 1 - loL
        loS <- floor(n * sumTrim) + 1
        hiS <- n + 1 - loS

#       keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
#       a fix from leonardo ivan almonacid cardenas, since rank() can return
#       non-integer values when there are a lot of ties

        keep <- (rank(logR)>=loL & rank(logR)<=hiL) & (rank(absE)>=loS & rank(absE)<=hiS)

        if(doWeighting)
                2^( sum(logR[keep]/v[keep], na.rm=TRUE) / sum(1/v[keep], na.rm=TRUE) )
        else
                2^( mean(logR[keep], na.rm=TRUE) )
}

