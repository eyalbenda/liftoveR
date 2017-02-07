# liftoveR:  liftover regions from one genome build to another by realignment

## Outline

The package *liftoveR* is a fast and versatile tools for moving between one set of genome coordinates to another based on sequence identity. This is done by realigning the sequence around the coordinates in the source genome to the target genome. The package incorporates multiple tools for alignment, including Bowtie (useful for stringently aligning closely related sequences), Rsubread (useful for aligning related sequences while allowing gaps) and blastn (mapping across species, or highly divergent sequences).
See the associated vignette for use instructions and use cases.

## Installation

The package depends on R packages that are available through Bioconductor:

* Biostrings
* GenomicRanges
* IRanges
* Rsamtools
* QuasR
*optional*: 
* Rsubread

To install these dependencies, run:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings","GenomicRanges","IRanges","Rsamtools","QuasR"))
```

To align using blastn, the package currently only support version *blast 2.2.31+*, available [here](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/). This is due to the SAM format output of blast, which changes with each version. If that settles down, I'll update to support their final version.


After that, the easiest way to install the package is if using the *devtools* package. If you want to quickly install the package, this should get it done.

```{r}
devtools::install_github("eyalbenda/liftoveR",args="--no-multiarch",build_vignette = FALSE)
```
*Note: Due to requirment of the package used for alignment, on windows only 64-bit versions of R are supported.*

The vignette include instructions to run an example dataset, using all three aligners. Since it includes performing alignments which are time consuming, building the *vignette* takes some time (~5 minutes, possibly more), so decide if it's worth it. To build with the vignette, run this:

```{r}
<<<<<<< HEAD
devtools::install_github("eyalbenda/liftoveR",args="--no-multiarch",build_vignette = TRUE)
=======
<<<<<<< HEAD
devtools::install_github("eyalbenda/liftoveR",args="--no-multiarch",build_vignette = TRUE)
=======
devtools::install_github("eyalbenda/liftoveR",args="--no-multiarch")
>>>>>>> 8611093599652eaf1e47167ca5887ecea41e1bb3

>>>>>>> b19965e9825a6ca3c5f12b6afa4b9201a744c6f7
```
Then, view it using
```{r}
vignette("liftoveR")
```
If needed, you can first install the *devtools* package using

```{r}
install.packages("devtools")
```
