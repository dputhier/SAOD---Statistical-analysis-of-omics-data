## Specify the URL of the base for the course
dir.base <- 'http://pedagogix-tagc.univ-mrs.fr/courses/statistics_bioinformatics'

## Load the general configuration file
source(file.path(dir.base, 'R-files', 'config.R'))

setwd(dir.results)
## Define the location of data directory and file containing expression profiles
dir.denpBer <- file.path(dir.base, 'data', 'gene_expression','denboer_2009')
file.expr.table <- file.path(dir.denboer, 'GSE13425_Norm_Whole.txt')

## Load the expression profiles (one row per gene, one column per sample).
##
## BEWARE, the whole file makes 20Mb.
## It can take a few minutes to be downloaded from the Web site and loaded in R.
expr.matrix <- read.table(file.expr.table, sep = "\t", head = T, row = 1)
print(dim(expr.matrix))
## Should give this: 22283   190
## Load the table describing each sample
## (one row per sample, one column per description field).
pheno <- read.table(file.path(dir.denboer, 'phenoData_GSE13425.tab'), sep='\t', head=TRUE, row=1)
dim(pheno)
## [1] 190   4

## We can check the content of the phenoData table by
## consulting its column names.
names(pheno)
## Print the values of the field Sample.title of pheno
print(pheno$Sample.title)
## Print the number of samples per cancer type
print(data.frame("n"=sort(table(pheno$Sample.title),decreasing=T)))
## Define an abbreviated name for each canceer subtype
## (will be used later visualization on the plots)
group.abbrev <- c(
      'BCR-ABL + hyperdiploidy'='Bch',
      'BCR-ABL'='Bc',
      'E2A-rearranged (E)'='BE',
      'E2A-rearranged (E-sub)'='BEs',
      'E2A-rearranged (EP)'='BEp',
      'MLL'='BM',
      'T-ALL'='T',
      'TEL-AML1 + hyperdiploidy'='Bth',
      'TEL-AML1'='Bt',
      'hyperdiploid'='Bh',
      'pre-B ALL'='pB'
                   )
sample.subtypes <- as.vector(pheno$Sample.title)
sample.labels <- group.abbrev[sample.subtypes]
names(sample.labels) <- names(expr.matrix)


## Check the label for a random selection of 10 samples.
## Each run should give a different result
sample(sample.labels, size=10)
## Define group-specific colors
group.colors <- c(
                    'BCR-ABL + hyperdiploidy'='cyan',
                    'BCR-ABL'='black',
                    'E2A-rearranged (E)'='darkgray',
                    'E2A-rearranged (E-sub)'='green',
                    'E2A-rearranged (EP)'='orange',
                    'MLL'='#444400',
                    'T-ALL'='violet',
                    'TEL-AML1 + hyperdiploidy'='#000066',
                    'TEL-AML1'='darkgreen',
                    'hyperdiploid'='red',
                    'pre-B ALL'='blue'
                    )

## Assign group-specific colors to patients
sample.colors <- group.colors[as.vector(pheno$Sample.title)]
names(sample.colors) <- names(expr.matrix)
sample(sample.colors,size=20)
## Plot the expression profiles of two arbitrarily selected genes
g1 <- 236
g2 <- 1213
x <- as.vector(as.matrix(expr.matrix[g1,]))
y <- as.vector(as.matrix(expr.matrix[g2,]))
dev.new(width=8,height=8)
plot(x,y,
       col=sample.colors,
       type="n",
       panel.first=grid(col="black"),
         main="Den Boer (2009), two selected genes",
         xlab=paste("gene", g1), ylab=paste("gene", g2))
text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend("topright",col=group.colors,
         legend=names(group.colors),pch=1,cex=0.6,bg="white",bty="o")
dev.copy2pdf(file="two_arbitrary_genes_plot.pdf", width=8, height=8)
## Compute sample-wise variance
var.per.sample <- apply(expr.matrix, 2, var)
head(var.per.sample)

## Inspect the distribution of sample-wise variance
dev.new(width=7,height=5)
hist(var.per.sample, breaks=20)

## Store the figure in a pdf file
dev.copy2pdf(file=file.path(dir.figures, "denboer2009_variance_per_sample.pdf"))
## Compute gene-wise variance
var.per.gene <- apply(expr.matrix, 1, var)

## Inspect the distribution of gene-wise variance
dev.new(width=7, height=5)
hist(var.per.gene, breaks=100)
dev.copy2pdf(file=file.path(dir.figures, "denboer2009_variance_per_gene.pdf"), width=7,height=5)
## Sort genes per decreasing variance
genes.by.decr.var <- sort(var.per.gene,decreasing=TRUE)

## Print the 5 genes with highest variance
head(genes.by.decr.var)
## Print the 5 genes with lowest variance
tail(genes.by.decr.var)
## Select the 30 top-ranking genes in the list sorted by variance.
## This list of genes will be used below as training variables for
## supervised classification.
top.nb <- 20 ## This number can be changed for testing
genes.selected.by.var <- names(genes.by.decr.var[1:top.nb])

## Check the names of the first selected genes
head(genes.selected.by.var, n=20)
## Create a data frame to store gene values and ranks
## for different selection criteria.
gene.ranks <- data.frame(var=var.per.gene)
head(gene.ranks)
## Beware, we rank according to minus variance, because
## we want to associate the lowest ranks to the highest variances
gene.ranks$var.rank <- rank(-gene.ranks$var, ties.method='random')
head(gene.ranks)
## Check the rank of the 5 genes with highest and lowest variance, resp.
gene.ranks[names(genes.by.decr.var[1:5]),]
## Print the bottom 5 genes of the variance-sorted table
gene.ranks[names(tail(genes.by.decr.var)),]
## Plot the expression profiles of the two genes with highest variance
maxvar.g1 <- names(genes.by.decr.var[1])
print(maxvar.g1)
(maxvar.g2 <- names(genes.by.decr.var[2]))
x <- as.vector(as.matrix(expr.matrix[maxvar.g1,]))
y <- as.vector(as.matrix(expr.matrix[maxvar.g2,]))
dev.new(width=8, height=8)
plot(x,y,
       col=sample.colors,
       type='n',
       panel.first=grid(col='black'),
       main="2 genes with the highest variance",
       xlab=paste('gene', maxvar.g1),
       ylab=paste('gene', maxvar.g2))
text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend('topright',col=group.colors,
         legend=names(group.colors),pch=1,cex=0.7,bg='white',bty='o')
dev.copy2pdf(file="denboer2009_2_maxvar_genes.pdf", width=8, height=8)
## Load a utility to apply Welch test on each row of a table
source(file.path(dir.util, "util_student_test_multi.R"))

## Define a vector indicating whether each sample
## belongs to the subtype of interest (e.g. "pB") or not.
current.group <- "pB"
one.group.vs.others<- sample.labels
one.group.vs.others[sample.labels != current.group] <- "other"
print(table(one.group.vs.others))
## Test the mean equality between pB subtype and all other subtypes
welch.one.group.vs.others <- t.test.multi(expr.matrix, one.group.vs.others)

## Store the volcano plot produced by the function t.test.multi()
dev.copy2pdf(file=file.path(dir.figures, paste(sep="", "denboer2009_welch_", current.group,"_vs_others_volcano.pdf")))

names(welch.one.group.vs.others)

## Update the gene rank table
test.name <- paste(current.group, '.vs.others.sig', sep='')
gene.ranks[,test.name] <- welch.one.group.vs.others$sig ## Add a column with significances
gene.ranks[,paste(test.name, ".rank", sep="")] <- rank(-welch.one.group.vs.others$sig, , ties.method='random') ## Add a column with significance ranks

## Apply the Welch test for the 3 other majority groups
for (current.group in c("Bh", "Bt", "T")) {
    one.group.vs.others<- sample.labels
    one.group.vs.others[sample.labels != current.group] <- "other"

    ## Test the mean equality between pB subtype and all other subtypes
    welch.one.group.vs.others <- t.test.multi(expr.matrix, one.group.vs.others)

    ## Store the volcano plot
    dev.copy2pdf(file=file.path(dir.figures, paste(sep="", "denboer2009_welch_", current.group,"_vs_others_volcano.pdf")))

    ## Update the gene rank table
    test.name <- paste(current.group, '.vs.others.sig', sep='')
    gene.ranks[,test.name] <- welch.one.group.vs.others$sig
    gene.ranks[,paste(test.name, ".rank", sep="")] <- rank(-welch.one.group.vs.others$sig, , ties.method='random')
}

## Check the resulting gene table
head(gene.ranks)
head(gene.ranks[order(gene.ranks$T.vs.others.sig.rank),])
tail(gene.ranks[order(gene.ranks$T.vs.others.sig.rank),])

## Store the gene rank table in a text file (tab-separated columns)
write.table(gene.ranks, file=file.path(dir.results, 'DenBoer_gene_ranks.tab'), sep='\t', quote=F, col.names=NA)
## Plot variance against significance of the Welch test for the 4 different groups
dev.new(width=8, height=8)
plot(gene.ranks[,c("var", "pB.vs.others.sig", "Bh.vs.others.sig", "Bt.vs.others.sig", "T.vs.others.sig")], col="grey")

## Plot ranks of variance against significance of the Welch test for the 4 different groups
dev.new(width=8, height=8)
plot(gene.ranks[,c("var.rank", "pB.vs.others.sig.rank", "Bh.vs.others.sig.rank", "Bt.vs.others.sig.rank", "T.vs.others.sig.rank")], col="grey")
dev.
## load the stats library to use the princomp() and prcomp() function
library(stats)

## Perform the PCA transformation
expr.prcomp <- prcomp(t(expr.matrix),cor=TRUE)

## Analyze the content of the prcomp result:
## the result of the method prcomp() is an object
## belonging to the class "prcomp"
class(expr.prcomp)
## Get the field names of the prcomp objects
names(expr.prcomp)

## Get the attributes of the prcomp objects
attributes(expr.prcomp)
dev.new(width=7, height=5)
plot(expr.prcomp, main='Den Boer (2009), Variance  per component', xlab='Component')
## Get the standard deviation and variance per principal component
sd.per.pc <- expr.prcomp$sdev
var.per.pc <- sd.per.pc^2

## Display the percentage of total variance explained by each
sd.per.pc.percent <- sd.per.pc/sum(sd.per.pc)
var.per.pc.percent <- var.per.pc/sum(var.per.pc)
barplot(var.per.pc.percent[1:10], main='Den Boer (2009), Percent of variance  per component', xlab='Component', ylab='Percent variance', col='#BBDDFF')
dev.copy2pdf(file="pca_variances.pdf", width=7, height=5)
dev.new(width=8, height=8)
biplot(expr.prcomp,var.axes=FALSE,
         panel.first=grid(col='black'),
         main=paste('PCA; Den Boer (2009); ',
           ncol(expr.matrix), 'samples *', nrow(expr.matrix), 'genes', sep=' '),
         xlab='First component', ylab='Second component')
## Plot components PC1 and PC2
dev.new(width=8, height=8)
plot(expr.prcomp$x[,1:2],
       col=sample.colors,
       type='n',
       panel.first=grid(col='black'),
         main=paste('PCA; Den Boer (2009); ',
           ncol(expr.matrix), 'samples *', nrow(expr.matrix), 'genes', sep=' '),
         xlab='PC1', ylab='PC2')
text(expr.prcomp$x[,1:2],labels=sample.labels,col=sample.colors,pch=0.5)
legend('bottomleft',col=group.colors,
         legend=names(group.colors),pch=1,cex=0.7,bg='white',bty='o')
dev.copy2pdf(file="PC1_vs_PC2.pdf", width=8, height=8)
## Plot components PC2 and PC3
plot(expr.prcomp$x[,2:3],
       col=sample.colors,
       type='n',
       panel.first=grid(col='black'),
         main=paste('PCA; Den Boer (2009); ',
           ncol(expr.matrix), 'samples *', nrow(expr.matrix), 'genes', sep=' '),
         xlab='PC2', ylab='PC3')
text(expr.prcomp$x[,2:3],labels=sample.labels,col=sample.colors,pch=0.5)
legend('bottomleft',col=group.colors,
         legend=names(group.colors),pch=1,cex=0.7,bg='white',bty='o')
dev.copy2pdf(file="PC2_vs_PC3.pdf", width=8, height=8)
## Load the library containing the linear discriminant analysis function
library(MASS)

## Redefined the variables (they may have changed during the previous processing)
current.group <- "pB"
one.group.vs.others<- sample.labels
one.group.vs.others[sample.labels != current.group] <- "other"
print(table(one.group.vs.others))

## Test the mean equality between pB subtype and all other subtypes
welch.pB.vs.others <- t.test.multi(expr.matrix, one.group.vs.others)

## Print the names of the 20 top-ranking genes.
## Note that these are not sorted !
print(rownames(welch.pB.vs.others[welch.pB.vs.others$rank <= 20,]))

## Sort gene names by Wekch sig
welch.pB.vs.others.sorted <- welch.pB.vs.others[order(welch.pB.vs.others$sig, decreasing=TRUE),]
sorted.names <- rownames(welch.pB.vs.others.sorted)
print(sorted.names[1:20])

## Print a plot with the two top-ranking genes with the pB versus other Welch test
g1 <- sorted.names[1]
g2 <- sorted.names[2]
x <- as.vector(as.matrix(expr.matrix[g1,]))
y <- as.vector(as.matrix(expr.matrix[g2,]))
     dev.new(width=8,height=8)
     plot(x,y,
     col=sample.colors,
     type="n",
     panel.first=grid(col="black"),
     main="Den Boer (2009), two selected genes",
     xlab=paste("gene", g1), ylab=paste("gene", g2))
text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend("topright",col=group.colors,
         legend=names(group.colors),pch=1,cex=0.6,bg="white",bty="o")

## Print a plot with the two top-ranking genes with the pB versus other Welch test
g3 <- sorted.names[3]
g4 <- sorted.names[4]
x <- as.vector(as.matrix(expr.matrix[g3,]))
y <- as.vector(as.matrix(expr.matrix[g4,]))
     dev.new(width=8,height=8)
     plot(x,y,
     col=sample.colors,
     type="n",
     panel.first=grid(col="black"),
     main="Den Boer (2009), two selected genes",
     xlab=paste("gene", g3), ylab=paste("gene", g4))
text(x, y,labels=sample.labels,col=sample.colors,pch=0.5)
legend("topright",col=group.colors,
         legend=names(group.colors),pch=1,cex=0.6,bg="white",bty="o")


## Select the 20 top-ranking genes sorted by decreasing variance
top.variables <- 4
selected.genes <- sorted.names[1:top.variables]


## Train the classifier
lda.classifier <- lda(t(expr.matrix[selected.genes,]),one.group.vs.others,CV=FALSE)
## Use the MASS:lda() function with the cross-validation option
lda.loo <- lda(t(expr.matrix[selected.genes,]),one.group.vs.others,CV=TRUE)

## Collect the LOO prediction result in a vector
loo.predicted.class <- as.vector(lda.loo$class)
print(loo.predicted.class)
table(loo.predicted.class)

## Build a contingency table of known versus predicted class
(lda.loo.xtab1 <- table(one.group.vs.others, loo.predicted.class))

## Check the assignation of all known groups to "pB" or "other"
(lda.loo.xtab2 <- table(sample.labels, loo.predicted.class))
## Display the contingency table as a heat map
image(lda.loo.xtab)
library(lattice)
levelplot(lda.loo.xtab)

## Compute the hit rate
hits <- sample.labels == loo.predicted.class
errors <- sample.labels != loo.predicted.class

## Compute the number of hits
## (we need to omit NA values because LDA fails to assign a group to some objects).
(nb.hits <- sum(na.omit(hits))) ## this should give 152
(nb.pred <- length(na.omit(hits))) ## This should give 187
(hit.rate <- nb.hits / nb.pred ) ## This should give 0.81
## Permute the training labels
sample.labels.perm <- as.vector(sample(sample.labels))

## Compare original training groups and permuted labels.
table(sample.labels, sample.labels.perm)
## Run LDA in cross-validation (LOO) mode with the permuted labels
lda.loo.labels.perm <- lda(t(expr.matrix[selected.genes,]),sample.labels.perm,CV=TRUE)

## Build a contingency table of known versus predicted class
loo.predicted.class.labels.perm <- as.vector(lda.loo.labels.perm$class)
lda.loo.labels.perm.xtab <- table(sample.labels.perm, loo.predicted.class.labels.perm)
print(lda.loo.labels.perm.xtab)
## Display the contingency table as a heat map
image(lda.loo.labels.perm.xtab)

## Levelplot uses a nice color palette
levelplot(lda.loo.labels.perm.xtab)

## Compute the number of hits
## (we need to omit NA values because LDA fails to assign a group to some objects).
hits.label.perm <- sample.labels.perm == loo.predicted.class.labels.perm
(nb.hits.label.perm <- sum(na.omit(hits.label.perm))) ## This gives 25 this time, but should give different results at each trial
(nb.pred.label.perm <- length(na.omit(hits.label.perm))) ## This should give 187
(hit.rate.label.perm <- nb.hits.label.perm / nb.pred.label.perm ) ## This gives 0.13 this time, should give different results at each trial
