## ----results="hide",warning=F,message=F----------------------------------
require(R.utils)
download.file("ftp://ftp.ensemblgenomes.org/pub/fungi/release-34/fasta/aspergillus_flavus/dna/Aspergillus_flavus.JCVI-afl1-v2.0.dna.toplevel.fa.gz","aspFlavusNew.fa.gz")
download.file("ftp://ftp.ensemblgenomes.org/pub/fungi/release-3/fasta/aspergillus_flavus/dna/Aspergillus_flavus.CADRE.55.dna.toplevel.fa.gz","aspFlavusOld.fa.gz")
gunzip("aspFlavusOld.fa.gz")
gunzip("aspFlavusNew.fa.gz")

## ----results="hide",warning=F,message=F----------------------------------
download.file("ftp://ftp.ensemblgenomes.org/pub/fungi/release-3/gtf/aspergillus_flavus/Aspergillus_flavus.CADRE.55.gtf.gz", "aspFlavusOld.gtf.gz")
download.file("ftp://ftp.ensemblgenomes.org/pub/fungi/release-34/gtf/aspergillus_flavus/Aspergillus_flavus.JCVI-afl1-v2.0.34.gtf.gz","aspFlavusNew.gtf.gz")
gunzip("aspFlavusOld.gtf.gz")
gunzip("aspFlavusNew.gtf.gz")

## ----results="hide",warning=F,message=F----------------------------------
require(stringr)
genesOldCoords = read.delim("aspFlavusOld.gtf",comment.char =  "#",header = F,as.is = T)
genesOldCoords = genesOldCoords[genesOldCoords[,3]=="exon",]

## ------------------------------------------------------------------------
head(genesOldCoords)

## ----results="hide",warning=F,message=F----------------------------------
geneNamesOld = str_replace_all(str_extract(genesOldCoords[,9],"gene_name+.[A-Z,a-z,0-9,_]*"),"gene_name ","")
exonNumberOld = str_replace_all(str_extract(genesOldCoords[,9],"exon_number+.[0-9]*"),"exon_number ","")
exonNameAndNumberOld = paste(geneNamesOld,exonNumberOld,sep="_")

## ------------------------------------------------------------------------
genesNewCoords = read.delim("aspFlavusNew.gtf",comment.char =  "#",header = F,as.is = T)
genesNewCoords = genesNewCoords[genesNewCoords[,3]=="exon",]
geneNamesNew = str_replace_all(str_extract(genesNewCoords[,9],"gene_id+.[A-Z,a-z,0-9,_]*"),"gene_id ","")
exonNumberNew = str_replace_all(str_extract(genesNewCoords[,9],"exon_number+.[0-9]*"),"exon_number ","")
exonNameAndNumberNew = paste(geneNamesNew,exonNumberNew,sep="_")

## ----results="hide",warning=F,message=F----------------------------------
exonsInBoth = intersect(exonNameAndNumberNew,exonNameAndNumberOld)
genesOldCoordsFinal = genesOldCoords[match(exonsInBoth,exonNameAndNumberOld),]
genesNewCoordsFinal = genesNewCoords[match(exonsInBoth,exonNameAndNumberNew),]

## ----message=F,warning=F-------------------------------------------------
require(liftoveR)
genesOldInNew = liftover(chrom = genesOldCoordsFinal[,1],
         start = genesOldCoordsFinal[,4],
         end = genesOldCoordsFinal[,5],
         originalBuild = "aspFlavusOld.fa",
         newBuild = "aspFlavusNew.fa",
         tmpdir = getwd(),
         maxMismatches = 5,
         aligner="Rbowtie")

## ------------------------------------------------------------------------
head(genesOldInNew)

## ------------------------------------------------------------------------
sum(is.na(genesOldInNew$chrom))

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom==genesNewCoords$V1 &genesOldInNew$start==genesNewCoords$V4 & genesOldInNew$end==genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom!=genesNewCoords$V1 |genesOldInNew$start!=genesNewCoords$V4 & genesOldInNew$end!=genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ----results="hide",warning=F,message=F----------------------------------
genesOldInNew = liftover(chrom = genesOldCoordsFinal[,1],
         start = genesOldCoordsFinal[,4],
         end = genesOldCoordsFinal[,5],
         originalBuild = "aspFlavusOld.fa",
         newBuild = "aspFlavusNew.fa",
         tmpdir=getwd(),
         maxMismatches = 5,
         aligner="Rsubread")

## ------------------------------------------------------------------------
head(genesOldInNew)

## ------------------------------------------------------------------------
sum(is.na(genesOldInNew$chrom))

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom==genesNewCoords$V1 &genesOldInNew$start==genesNewCoords$V4 & genesOldInNew$end==genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom!=genesNewCoords$V1 |genesOldInNew$start!=genesNewCoords$V4 & genesOldInNew$end!=genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ----results="hide",warning=F,message=F----------------------------------
genesOldInNew = liftover(chrom = genesOldCoordsFinal[,1],
         start = genesOldCoordsFinal[,4],
         end = genesOldCoordsFinal[,5],
         originalBuild = "aspFlavusOld.fa",
         newBuild = "aspFlavusNew.fa",
         tmpdir=getwd(),
         lengthSides = 100,
         maxMismatches = 30,
         aligner="Rsubread")

## ------------------------------------------------------------------------
sum(is.na(genesOldInNew$chrom))

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom==genesNewCoords$V1 &genesOldInNew$start==genesNewCoords$V4 & genesOldInNew$end==genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom!=genesNewCoords$V1 |genesOldInNew$start!=genesNewCoords$V4 & genesOldInNew$end!=genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100


## ------------------------------------------------------------------------
registerBlast("/usr/bin/blastn","/usr/bin/makeblastdb","/usr/bin/samtools")

## ------------------------------------------------------------------------
genesOldInNew = liftover(chrom = genesOldCoordsFinal[,1],
         start = genesOldCoordsFinal[,4],
         end = genesOldCoordsFinal[,5],
         originalBuild = "aspFlavusOld.fa",
         newBuild = "aspFlavusNew.fa",
         tmpdir=getwd(),
         aligner="blastn")

## ------------------------------------------------------------------------
sum(is.na(genesOldInNew$chrom))

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom==genesNewCoords$V1 &genesOldInNew$start==genesNewCoords$V4 & genesOldInNew$end==genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ------------------------------------------------------------------------
sum(genesOldInNew$chrom!=genesNewCoords$V1 |genesOldInNew$start!=genesNewCoords$V4 & genesOldInNew$end!=genesNewCoords$V5,na.rm = T) / nrow(genesOldInNew) * 100

## ----results="hide",warning=F,message=F----------------------------------
removeDirectory("aspFlavusNew.fa.Rbowtie",recursive = T)
tmpdir = dir(pattern="Rtmp")
removeDirectory(tmpdir,recursive = T)
file.remove(dir(pattern="aspFlavus"))
file.remove(dir(pattern="QuasR"))
file.remove(dir(pattern="liftOverInput"))
file.remove(dir(pattern="blast"))

