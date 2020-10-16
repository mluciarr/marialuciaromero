# marialuciaromero
**ANALYSIS OF MARMOSET mRNA SEQUENCES**
===================================
*by* ROMERO, Maria Lucia


# Introduction

We are going to reuse an existing report, describing a simple analysis protocol, to analyze a new dataset.

You can submit your worked exercise to the following [UPF Virtual Campus link](https://aulaglobal.upf.edu/mod/assign/view.php?id=39525). Please take a look there for further instructions about the submission.

## Objectives

* This exercise show many of the strengths of the Unix command-line to process genomic data.
* We will also practice reporting such commands and the results using a text file and the `MarkDown` syntax.

First of all we need to create a project folder, of course, and set it as the current working directory...

```{.sh}
mkdir exercise
cd exercise
```

# Discussion 

We can observe that in the [Table 1](#table1) that the median value in the nucleotide length analysis is slightly more significant than the mean because the minimum value and the maximum are strongly diferente. On the other part, the values of media and mean in the GC content are quite similar due to having a normal distribution.

The plot named [Figure 2](#fig2) shows that the frequency of GC percentage has a normal distribution of the data. Where there are more mRNA with a GC content around 54 percent. On the other hand, the plot of frequency of nucleotide length ([Fig. 1](#fig1)) shows us that the vast majority of the mRNA length is around 680 and 708 bp.

In [Figure 4](#fig4) we can analyse the relationship between nucleotide length and GC content. Where it has the GC content histogram at the top of the plot which below to the "x" axis  and the nucleotide length on the right of the graph which below to "y" axis, in order to see whether or not there is a relationship between the length of mRNA the the concentration of GC. We have observed that the relationship has not a normal distribution. However we can see that there is more concentration of the data between the 40 and 60 percent of GC and from 0 to 1000 of bp (base pare) in "y axis". We can not find a relationship between the length of the mRNA and the percentage of GC, due to most of the mRNA length from 41 to 708 have a wide range of contraction of GC. Therefore all  size of mRNA can have from 30 to 70 percent of GC.

As the trend-line of [Figure 4](#fig4) shows, there is not a significant correlation between the length of the nucleotides and the GC content in them due to that  the sizes of mRNA (around 700bp) remains content for different concentration of GC. However, we can notice that from  50 percentage of GC concentration there is a slight increase of length. This could means that the more concentration of GC the longer are the mRNA. Although, we cannot ensure that  increase is real  due to not be significant.

Finally, it would be interesting to make a deeper study about where they are located in chromosomes and their efficient transcription to analyze if mRNA processing is responsible for the high expression of GC-rich genes.


# Protocol

## The original data

We have to analyze sequence length distribution and GC content for the current mRNA sequences annotated on *Callithrix jacchus* genome (March, 2009). We have connected to the [UCSC genome browser download web site](http://hgdownload.soe.ucsc.edu/downloads.html), and followed the [Full data set link](http://hgdownload.soe.ucsc.edu/goldenPath/calJac3/bigZips/). We can see that there is a file named `mrna.fa.gz`, we can just copy the corresponding link to the following command:

```{.sh}
ftp://hgdownload.cse.ucsc.edu/goldenPath/calJac3/bigZips/mrna.fa.gz

```

We obtained a compressed file containing all the mRNA sequences annotated for this genome in `fasta` format. We can check it's content:

```{.sh}
zcat <  mrna.fa.gz | head -10
# > EF214838 1
# gaggcttgcaaagtgaatccttttgtcaaggacaggtggggcaacattcc
# cctggatgatgctgtacagttcaaccatctggaggtagtcaaactgcttc
# aagattaccaggactcctacacactctctgagactcaggctgaggcagca
# gctgaagctctgtccaaagagaacttagagagcatggtatgagcacaggt
# catggatagcctctgctcaagaaaaagcattagctggccacacatgcaat
# ccatacccaccaaaaatgctatggagaactatactgcttcagtcgggacc
# aagcagtcatttggtgacttggggtactgctttctatgagagtcaaaata
# ccccattccctcagcagacagagtacagagaagggctatctctggatagc
# agtagagctatccagagagaatgggcttcaaagtacagccaaatggcttg

```

Let's see how many mRNA sequences do we have. We can uncompress the file and keep going with the flat text file, which can take a lot of disk resources, or we can keep using the compressed file as in the previous command.

```{.sh}
zcat < mrna.fa.gz | egrep -c '^>' 
# 1596
```

In case that a program can deal with compressed files, we can take advantage of that feature. There is an `egrep` version that can read such files, if you were guessing is `zegrep` of course (as it happens with `cat` and `zcat`).

```{.sh}
zegrep -c '^>' mrna.fa.gz
# 1596
```

## Analysis of CDS sequences

To calculate the nucleotide length and the GC content, we can use on e of the programs that is provided within the [EMBOSS suite](http://emboss.sourceforge.net/): `infoseq`. You can get information about this tool [from this link](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/infoseq.html).

```{.sh}
infoseq -help
# you get info about this tool command-line options
# you can also use: man infoseq

zcat < mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -noheading -only -name -length -pgc | \
    head -5
# Display basic information about sequences
# EF214838       680    47.65  
# EF214839       680    50.59  
# EF214840       680    53.09  
# EF214841       680    51.47  
# EF214842       680    35.74  
```

In the above command we are just looking to the expected output, the following is doing the job and saving the output into `mrna.lengc.tbl` file.

```{.sh}
zcat < mrna.fa.gz | \
  infoseq -sequence fasta::stdin \
          -outfile mrna.lengc.tbl \
          -noheading -only -name -length -pgc
```

## Visualizing the analysis

By running `R` command, we enter in the `R` shell interpreter, which understands `R` commands

```{.sh}
R
# 
# R version 3.5.1 (2018-07-02) -- "Feather Spray"
# Copyright (C) 2018 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
```

Now, we must load the tabular data into a variable.

```{.r}
DATA <- read.table("mrna.lengc.tbl", header=FALSE);

# just checking the data structure
head(DATA,4);
#         V1  V2    V3
# 1 EF214838 680 47.65
# 2 EF214839 680 50.59
# 3 EF214840 680 53.09
# 4 EF214841 680 51.47
```

We can rename the table columns, so they are more meaningful:

```{.r}
#colnames(DATA) <- c("ID","NUClen","GCpct");
head(DATA,4);
#         ID NUClen GCpct
# 1 EF214838    680 47.65
# 2 EF214839    680 50.59
# 3 EF214840    680 53.09
# 4 EF214841    680 51.47
```

Let's calculate some stats on the dataset.

<a name="table1"> *Table 1*: Statistics analysis of the data </a>

```{.r}
summary(DATA)
#         ID           NUClen            GCpct
#  AB038384:   1   Min.   :  41.0   Min.   :28.07  
#  AB038385:   1   1st Qu.: 450.8   1st Qu.:46.52  
#  AB038386:   1   Median : 680.0   Median :54.55  
#  AB046546:   1   Mean   : 708.1   Mean   :53.11  
#  AB046547:   1   3rd Qu.: 789.8   3rd Qu.:60.01  
#  AB046548:   1   Max.   :9396.0   Max.   :73.82  
#  (Other) :1590                             ```
```

Now, let's make an histogram:

```{r.}
hist(DATA$NUClen);
hist(DATA$GCpct);
```


<a name="fig1">![Showing the frequency of GC content](/Users/LuciaR/Documents/Bioinformatics_4_health_science/ALG/markdown/exercise/RplotDATA$NUClen.png "Showing the frequency of Nucleotide length")</a>

_Figure 1_: Showing the frequency of nucleotide length.

<a name="fig2">![Showing the frequency of GC content](/Users/LuciaR/Documents/Bioinformatics_4_health_science/ALG/markdown/exercise/RplotDATA$GCpct.png "Showing the frequency of GC content")</a>


_Figure 2_: Showing the frequency of GC content.

Or to compare both measures and save into a PNG image...

```{.r}
png(file="plot.png");
plot(DATA$NUClen ~ DATA$GCpct);
dev.off();
```

<a name="fig3"> ![Showing GC content versus sequence length](/Users/LuciaR/Documents/Bioinformatics_4_health_science/ALG/markdown/exercise/plot.png "Showing GC content versus sequence length")</a>


_Figure 3_: Showing GC content versus sequence length


Just take care that this is not a single-variable distribution but a two axis scatter-plot, each measure is a continuous random variable that may fit (or not) a normal distribution. The real data distribution is shown on the histograms, we can merge then those three plots...

```{.r}
png(file="plot2.png");
def.par <- par();
# preparing a layout grid where to combine different plots
nf <- layout(matrix(c(2,0,1,3), # matrix contents is plot order
                    2, 2, byrow=TRUE),
             c(3,1), c(1,3), # cell relative sizes
             TRUE);
# computing data distribution
xhist <- hist(DATA$GCpct, plot=FALSE);
yhist <- hist(DATA$NUClen, plot=FALSE);
datamax <- max(xhist$counts,
               xhist$counts);
#
# drawing the main plot
par(mar=c(5,5,1,1)); 
plot(DATA$NUClen ~ DATA$GCpct,
      main="",
      xlab="GC%",
      ylab="Sequence length (bp)",
      col="seagreen2");
lines(lowess(DATA$NUClen ~ DATA$GCpct),
      col="red", lty="dotted");
mtext(paste("n=",nrow(DATA),sep=""),
      side=3, line=-1);
# drawing x-axis histogram
par(mar=c(0,5,1,1));
barplot(xhist$counts, ylim=c(0, datamax),
        axes=FALSE, space=0, 
        col="indianred3", main="C.elegans mRNAs")
# drawing y-axis histogram
par(mar=c(5,0,1,1));
barplot(yhist$counts, xlim=c(0, datamax),
        axes=FALSE, horiz=TRUE, space=0, 
        col="indianred3", main="")
#
par(def.par); # reseting graphical parameters
dev.off();
```

<a name="fig4">![Showing GC content versus sequence length and marginal distributions\label{fig:GC&NUClen}](/Users/LuciaR/Documents/Bioinformatics_4_health_science/ALG/markdown/exercise/plot2.png "Showing GC content versus sequence length and marginal distributions")</a>

*Figure 4:* Showing GC content versus sequence length and marginal distributions


# CODA

## Software

We have used the following versions:

```{.sh}
uname -a
# Darwin MacBook-Pro-de-Lucia.local 19.6.0 Darwin Kernel Version 19.6.0: Mon Aug 31 22:12:52 
PDT 2020; root:xnu-6153.141.2~1/RELEASE_X86_64 x86_64

wget --version
# GNU Wget 1.20.3 hecho en darwin19.0.0.

infoseq -version
# EMBOSS:6.6.0.0

R --version
# R version 3.5.1 (2018-07-02)-- "Eggshell Igloo"
# Copyright (C) 2018 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin15.6.0 (64-bit)

```

## About this document

This document can be compiled into a PDF with the following command:

```{.sh}
pandoc -f gfm -t latex                  \
       --pdf-engine=xelatex             \
       --variable papersize:a4paper     \
       --variable geometry:margin=1.5cm \
       --variable fontsize=11pt         \
       --highlight-style pygments       \
       -o ROMERO_MARIALUCIA_ALG2021_genomicanalysis.pdf    \
          ROMERO_MARIALUCIA_ALG2021_genomicanalysis.md

```

You can get further information from the following links about the [MarkDown syntax](http://daringfireball.net/projects/markdown/syntax#link) and the [GitHub flavour](https://help.github.com/articles/github-flavored-markdown/). The code above requires `pandoc` and some `LaTeX` packages installed in your system, you can use `synaptic`, `apt-get` or `aptitude` to retrieve and install those tools from linux repositories. **Check out** that in the original text document we use three *back-tildes* to delimit the code blocks, to distinguish from the rest of the report.
