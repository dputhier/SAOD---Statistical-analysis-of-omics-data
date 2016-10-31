
## ------------------------------------------------------------------------
g <- 75 ## Number of submitted genes
k <- 59 ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
m <- 611 ## Number of "marked" elements, i.e. genes associated to this biological process
N <- 13588 ## Total number of genes with some annotation in GOTERM_BP_FAT.  
n <- N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
x <- 19 ## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process


## ------------------------------------------------------------------------

## Percent of the gene selection involved in "cell cycle". This
## corresponds to the "%" column returned by DAVID.
(percent.selection <- x / g * 100)

## Percent, among submitted genes with at least one annotation, which
## are involved in "cell cycle".
(percent.selection <- x / k * 100) 

## Percent of the biological process covered by the gene set.
(percent.process <- x / m * 100) 

## Random expectation
(marked.proportion <- m / N)
(exp.x <- k * marked.proportion)

## Fold enrichment, as computed by David
(fold.enrichment <-  (x / k ) / (m / N)) 
    


## ------------------------------------------------------------------------
## Define the range of possible values for x.
## The number of marked elements in the selection can neither be 
## higher than the size of the selection, nor higher than the 
## number of marked elements.
x.range <- 0:min(k,m) 

help(dhyper) ## Always read the documentation of a function before using it. 

## Compute the distribution of density P(X=x)
dens <- dhyper(x=x.range, m=m, n=n, k=k)

## Plot the distributon of hypergeometric densities
plot (x.range, dens, type="h", lwd=2, col="blue", main="Hypergeometric density", xlab="x = marked elements in the selection", ylab="density = P(X=x)", ylim=c(0, max(dens)*1.25))

## Draw an arrow indicating the expected value
arrows(exp.x, max(dens)*1.15, exp.x, max(dens)*1.05, lwd=2, col='darkgreen', angle=30, length=0.1, code=2)
text(exp.x, max(dens)*1.20, labels=paste("exp=", round(exp.x, digits=2)), col="darkgreen", font=2)

## Draw an arrow indicating the observed value
arrows(x, max(dens)*1.15, x, max(dens)*1.05, lwd=2, col='red', angle=30, length=0.1, code=2)
text(x, max(dens)*1.20, labels=paste("x=", x), col="red", font=2)
  


## ------------------------------------------------------------------------

## Plot the p-value distribution with logarithmic axes
plot (x.range, dhyper(x=x.range, m=m, n=n, k=k), type="l", lwd=2, col="blue", main="Hypergeometric density (log Y scale)", xlab="x = marked elements in the selection", ylab="density = P(X=x)", log="y", panel.first=grid())

## Arrow indicating expected value
arrows(exp.x, max(dens)*1e-10, exp.x, max(dens)*1e-2, lwd=2, col='darkgreen', angle=30, length=0.1, code=2)
text(exp.x, max(dens)*1e-15, labels=paste("exp=", round(exp.x, digits=2)), col="darkgreen", font=2)

## Arrow indicating observed value
arrows(x, dens[x+1]*1e-10, x, dens[x+1]*1e-2, lwd=2, col='red', angle=30, length=0.1, code=2)
text(x, dens[x+1]*1e-15, labels=paste("x=", x), col="red", font=2)




## ------------------------------------------------------------------------
print(dhyper(x=k, m=m, n=n, k=k))


## ------------------------------------------------------------------------
p.value <-  phyper(q=x -1, m=m, n=n, k=k, lower.tail=FALSE)


## ------------------------------------------------------------------------
      help(phyper)
    


## ------------------------------------------------------------------------
## Compute the distribution of P-values (for all possible values of x).
p.values <- phyper(q=x.range -1, m=m, n=n, k=k, lower.tail=FALSE)

## Plot the P-value distribution on linear scales
plot(x.range, p.values, type="l", lwd=2, col="violet", main="Hypergeometric P-value", xlab="x = marked elements in the selection", ylab="P-value = P(X>=x)", ylim=c(0, 1), panel.first=grid())

## Arrow indicating observed value
arrows(x, max(dens)*1.35, x, max(dens)*1.1, lwd=2, col='red', angle=30, length=0.1, code=2)
text(x, max(dens)*1.5, labels=paste("x=", x, "; p-val=", signif(digits=1, p.value), sep=""), col="red", font=2)

## We can plot the density below the P-value
lines(x.range, dens, type="h", col="blue", lwd=2)
legend("topright", legend=c("P-value", "density"), lwd=2, col=c("violet", "blue"), bg="white", bty="o")



## ------------------------------------------------------------------------
## Plot the P-value distribution with a logarithmic Y scale

plot(x.range, p.values,  type="l", lwd=2, col="violet", main="Hypergeometric P-value (log Y scale)", xlab="x = marked elements in the selection", ylab="P-value = P(X>=x); log scale", panel.first=grid(), log='y')

## Arrow indicating observed value
arrows(x, p.value, x, min(dens), lwd=2, col='red', angle=30, length=0.1, code=2)
text(x, min(dens), labels=paste("x=", x, sep=""), col="red", font=2, pos=4)

## Arrow indicating the P-value
arrows(x, p.value, 0, p.value, lwd=2, col='red', angle=30, length=0.1, code=2)
text(0, p.value*1e-5, labels=paste("p-val=", signif(digits=2, p.value), sep=""), col="red", font=2, pos=4)




## ------------------------------------------------------------------------

## Prepare a two-dimensional contingency table
contingency.table <- data.frame(matrix(nrow=2, ncol=2))
rownames(contingency.table) <- c("predicted.target", "non.predicted")
colnames(contingency.table) <- c("class.member", "non.member")

## Assign the values one by one to make sure we put them in the right
## place (this is not necessary, we could enter the 4 values in a
## single instruction).
contingency.table["predicted.target", "class.member"] <- x ## Number of marked genes in the selection
contingency.table["predicted.target", "non.member"] <- k - x ## Number of non-marked genes in the selection
contingency.table["non.predicted", "class.member"] <- m - x ## Number of marked genes outside of the selection
contingency.table["non.predicted", "non.member"] <- n - (k - x) ## Number of non-marked genes in the selection


print(contingency.table)


## Print marginal sums
(contingency.row.sum <- apply(contingency.table, 1, sum))
(contingency.col.sum <- apply(contingency.table, 2, sum))

## Create a contingency table with marginal sums
contingency.table.margins <- cbind(contingency.table, contingency.row.sum)
contingency.table.margins <- rbind(contingency.table.margins, apply(contingency.table.margins, 2, sum))
names(contingency.table.margins) <- c(names(contingency.table), "total")
rownames(contingency.table.margins) <- c(rownames(contingency.table), "total")
print(contingency.table.margins)


## Check the total
print(sum(contingency.table)) ## The value shoudl equal N, since every
                              ## possible gene must be assigned to one
                              ## cell of the contingency table.
print(N)

## Run Fisher's exact test
ftest.result <- fisher.test(x=contingency.table, alternative="greater")
print(ftest.result)
attributes(ftest.result) ## Display the list of attribute of the object returned by ftest
print (ftest.result$p.value) ## Print the P-value of the exact test

################################################################
## Compute expected values in the contingency table
exp.contingency.table <- contingency.table ## Quick and dirty way to obtain the same structure as the original contingency table
exp.contingency.table[ ] <- contingency.row.sum %*% t(contingency.col.sum) / N

## Add row and column sums (marginal sums)
exp.contingency.table <- cbind(exp.contingency.table, contingency.row.sum)
exp.contingency.table <- rbind(exp.contingency.table, apply(exp.contingency.table, 2, sum))
names(exp.contingency.table) <- c(names(contingency.table), "total")
rownames(exp.contingency.table) <- c(rownames(contingency.table), "total")
print(exp.contingency.table)
  


