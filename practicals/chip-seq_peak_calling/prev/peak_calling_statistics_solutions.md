Statistics for bioinformatics - Exercises - Peak calling statistics
===================================================================

* * * * *

Objectives
----------

During this practical, we will develop in R a simplified version of the
MACS peak calling algorithm, and test it on a relatively small dataset
(ChIP-seq data for the FNR transcription factor in the genome of the
Bacteria Escherichia coli).

Beware: the resulting script will certainly not be as good as existing
peak calling algorithms. The essential goal of this tutorial is to get a
deeper insight into the statistics behind peak calling, by implementing
ourselves a simplified version of the statistical test.

* * * * *

Data sets
---------

For each sample (FNR ChIP-seq and genomic input) we provide bed files
indicating the number of reads mapped in non-overlapping windows of
200bp, or 50bp, respectively (two files per sample).

  ------------------------------------ ------------------------------------
  FNR ChIP-seq sample (test)           FNR\_ChIP-seq\_Anaerobic\_A\_GSM1010
                                       219\_reads\_per\_200bp.bedg
                                       FNR\_ChIP-seq\_Anaerobic\_A\_GSM1010
                                       219\_reads\_per\_50bp.bedg

  Genomic input (control)              Escherichia\_coli\_K\_12\_MG1655\_in
                                       put\_reads\_per\_200bp.bedg
                                       Escherichia\_coli\_K\_12\_MG1655\_in
                                       put\_reads\_per\_50bp.bedg
  ------------------------------------ ------------------------------------

These files are in [bedGraph
format](http://genome.ucsc.edu/goldenPath/help/bedgraph.html).

**Remember:** the bed convention uses zero-based coordinates, with
semi-open intervals. Thus, the coordinates `0  50` correspond to

-   the semi-open interval `[0:50[` in zero-based coordinatesl
-   i.e. the closed interval `[0:49]` in zero-based coordinates;
-   i.e. the closed interval `[1:50]` in one-based (human
    understandable) coordinates.

### How was this dataset generated?

This section is optional. It explains the tricks I used to generate a
data set for this exercise.

I obtained the original datasets (mapped reads for the ChIP-seq and
input samples) by following **Morgane Thomas-Chollier's** tutorial
[Hands-on introduction to ChIP-seq
analysis](http://www.biologie.ens.fr/~mthomas/other/chip-seq-training/).

I then used bedtools to count the number of reads per bins, i.e.
non-overlapping windows of fixed size (200bp per window) covering the
full genome of Escherichia coli.

The protocol to count reads per bin is the following:

1.  Generate a bed file defining the bin limits (one row per bin)
2.  Convert the mapped reads from sam to bed file, sorted by chromosomal
    position
3.  Use bedtools to compute the intersection between each bin and the
    read file (i.e. count the reads falling onto each bin)

#### Generating windows of equal size along the reference genome

``` {.brush:bash;}
      bedtools makewindows -g Escherichia_coli_K_12_MG1655_genome.txt -w 200 \
      > Escherichia_coli_K_12_MG1655_windows_200bp.bed
    
```

The file genome.txt is expected to be a 2-columns file indicating the ID
and length of each chromosome. In the case of Escherichia coli, the file
contains a single line:

    gi|49175990|ref|NC_000913.2|    4639675

#### Converting bam to sorted bed

To compute the intersection between two sets of genomic regions
(bedtools intersect), bedtools requires two bed files sorted by
chromosome and by chromosomal position.

``` {.brush:bash;}
      ## We first need ton convert sam to bed format.
      ## For this I use an intermediate BAM format (I foudn a tool sam2bed but it does not seem to work)

      ## Convert mapped reads of FNR ChIP-seq library
      samtools view -bS SRR576933.sam | bedtools bamtobed | sort -n -k 2 > SRR576933_sorted.bed

      ## Convert mapped reads of control library
      samtools view -bS SRR576938.sam | bedtools bamtobed | sort -n -k 2 > SRR576938_sorted.bed
    
```

#### Counting the reads per bin

``` {.brush:bash;}
      ## Count the number of reads per window in the FNR ChIP-seq library
      bedtools intersect -a Escherichia_coli_K_12_MG1655_windows_200bp.bed \
      -b SRR576933_sorted.bed \
      -c -sorted \
      >  FNR_ChIP-seq_Anaerobic_A_GSM1010219_reads_per_200bp.bedg

      ## Count the number of reads per window in the control library
      bedtools intersect -a Escherichia_coli_K_12_MG1655_windows_200bp.bed \
      -b SRR576938_sorted.bed \
      -c -sorted \
      >  Escherichia_coli_K_12_MG1655_input_reads_per_200bp.bedg
    
```

* * * * *

Questions
---------

1.  Download the two bed files describing read maps.
2.  Load these two tables in R.
3.  Count the total reads for the FNR and input libraries, respectively.
4.  Normalize the input library in order to obtain the same sum as the
    test (FNR) library.
5.  Compare the distributions of counts per reads (normalized for the
    input) in the ChIP-seq and input samples.

6.  For each bin, compute the following statistics, and collect them in
    a table (one row per bin, one column per statistics):

    **Note**: each of these statistics can be computed with a single
    operation -- in R, you should avoid loops whenever possible.

    1.  Number of reads in the test (FNR ChIP-seq)
    2.  Number of reads in the input
    3.  Normalized number of reads in the input (norm.input)/
    4.  Reads per base in the input
    5.  Fold enrichment (ratio between test and normalized input)
    6.  Log-fold enrichment (log~10~ of the fold enrichment)
    7.  P-value of the test reads, using a Poisson distribution with
        local lambda estimate (explain your model).
    8.  P-value of the test reads, using a binomial distribution
        (explain your model).
    9.  P-value of the test reads, using a hypergeometric distribution
        (explain your model).

7.  Draw some dot plots to compare the different statistics (fold
    enrichment, log-fold, p-values with the different models).
8.  **Discuss the results.**

* * * * *

References
----------

1.  MACS method descripion: \
    Zhang, Y., Liu, T., Meyer, C.A., Eeckhoute, J., Johnson, D.S.,
    Bernstein, B.E., Nussbaum, C., Myers, R.M., Brown, M., Li, W., et
    al. (2008) Model-based analysis of ChIP-Seq (MACS). Genome Biol, 9,
    R137.

2.  To generate the data, we followed **Morgane Thomas-Chollier's**
    tutorial [Hands-on introduction to ChIP-seq
    analysis](http://www.biologie.ens.fr/~mthomas/other/chip-seq-training/),
    and then applied some tricks to count reads in fixed-width windows
    over the whole genome of *Escherichia coli*.

* * * * *
