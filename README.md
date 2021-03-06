# liftoveR:  liftover regions from one genome build to another by realignment

## Outline

The package *liftoveR* is a fast and versatile tools for moving between one set of genome coordinates to another based on sequence identity. This is done by realigning the sequence around the coordinates in the source genome to the target genome. The package incorporates multiple tools for alignment, including Bowtie (useful for stringently aligning closely related sequences), Rsubread (useful for aligning related sequences while allowing gaps) and blastn (mapping across species, or highly divergent sequences).


## Installation

The package depends on R packages that are available through Bioconductor:

* Biostrings
* GenomicRanges
* IRanges
* Rsamtools
* RBowtie
*optional*: 
* Rsubread

To install these dependencies, run:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings","GenomicRanges","IRanges","Rsamtools","RBowtie"))
```

To align using blastn, the package currently only support version *blast 2.2.31+*, available [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/). This is due to the SAM format output of blast, which changes with each version. If that settles down, we'll update to support their final version.


After that, the easiest way to install the package is if using the *devtools* package. If you want to quickly install the package, this should get it done.

```{r}
devtools::install_github("eyalbenda/liftoveR",args="--no-multiarch")
```
*Note: Due to requirment of the package used for alignment, on windows only 64-bit versions of R are supported.*

If needed, you can first install the *devtools* package using

```{r}
install.packages("devtools")
```
