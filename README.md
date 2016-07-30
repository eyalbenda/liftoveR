# liftoveR:  Lift over regions from one genome build to another by realignment

## Outline

Everyone who has worked with genomic data knows the frustration of working with different builds, and trying to fit everything together. To help with this, I wrote the package *liftoveR*. The package is highly versatile, and useful for any organism. All you need are the coordinates of your variants, the source genome, and the target you want to lift over onto. The package is made possible by the R community, and specifically packages for alignment of short reads from the *Bioconductor* repository (see next section for installation instructions).

## Installation

The package depends on R packages that are available through Bioconductor:

* Biostrings
* GenomicRanges
* IRanges
* Rsamtools
* QuasR

To install these dependencies, run:

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings","GenomicRanges","IRanges","Rsamtools","QuasR"))
```


After that, the easiest way to install the package is if using the *devtools* package.

```{r}
devtools::install_github("eyalbenda/liftoveR",args="--no-multiarch")
```
*Note: Due to requirment of the package used for alignment, on windows only 64-bit versions of R are supported.*


If needed, you can first install the *devtools* package using

```{r}
install.packages("devtools")
```

## Usage example and notes

All the work is done by the liftover function.
For example, if we have the following variants:

| chromosome    | start |  end |
| ------------- | ------| -----|
| chrI          | 85672 | 85672|
| chrII         | 6533  |6533 |
| chrIII        | 84430 |84430 |
| chrIV         | 4779  |4779  |

  We also have two fasta files, the original build to which these variants match:

  ```{r}
originalGenome = "originalGenome.fa"
```

And the build we want to lift over to:


  ```{r}
newGenome = "newGenome.fa"
```

Lift over can be done by a simple one liner:
  ```{r}
newCoordinates = liftover(chromosome,start,end,originalGenome,newGenome)
```
If you don't supply feature names, the package will create them by itself in the format *chromosome_start_end*.

Also, fasta index files (".bai files") will be generated automatically if missing,  and a Bowtie index for alignment to the new build is generated in a subdirectory called RBowtie. These files are currently not deleted by the package!

**Note!** Currently, due to the features of the package QuasR, lift over generates a log file in the R working directory related to the realignment.
