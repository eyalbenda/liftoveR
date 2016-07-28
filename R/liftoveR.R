#' @import QuasR
#' @import Rsamtools
#' @import Biostrings
#' @import GenomicRanges
#' @import IRanges
#' @importFrom utils write.table
indexFasta = function(fa)
{
  sprintf("Building fasta index for %s",fa)
  if(!file.exists(sprintf("%s.fai",fa)))
    indexFa(fa)
}

checkFasta = function(fa)
{
  firstLine = read.table(fa,nrows=1,stringsAsFactors = F)
  if(!substr(firstLine,1,1)[1] == ">")
    return(F)
  T
}

#' Lift over coordinates of variants from one genome build to another, by realignment
#' @param chrom For each variant, the chromosome on which it resides
#' @param start,end Physical coordinates of variants (in 1-based indexing. see \url{http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms})
#' @param originalBuild Path to the fasta sequence of the build from which the variants originate
#' @param newBuild Path to the fasta sequence of the build to which the variants will be lifted over
#' @param varnames Variant names. If missing, variants will be given names based as chromosome_position in original build
#' @param lengthSides number of bases which will be extracted from each side of the variant
#' @return A data frame containing the variants with their positions in the new build, based on the following scheme:
#' \describe{
#' \item{name}{Name of the variant}
#' \item{chrom_orig, start_orig, end_orig}{The original position supplied to function.}
#' \item{chrom, start, end}{Coordinates of alignment in the new build}
#' }
#' \strong{Note:} Currently, a log file is saved to the current working directory. This is a feature of the QuasR package used for mapping, and it is currently not possible to disable it.
#' @export
liftover = function(chrom,start,end = NULL,originalBuild,newBuild,varnames=NULL,lengthSides=25,tmpdir = NULL)
{
  if(!is.numeric(lengthSides)|length(lengthSides)!=1|lengthSides<1|round(lengthSides)!=lengthSides) stop("Error: lengthSides must be a single positive round number")
  if(is.null(end)) end = start
  if(!is.numeric(start)|!is.numeric(end)) stop("Error: Start and/or End positions supplied are not numeric vectors!")
  if(!all(end>=start)) stop("Error: all End positions must be at >= Start positions")
  if(!file.exists(originalBuild)|!file.exists(newBuild)) stop("Error: originalBuild and/or newBuild are not paths to existing files")
  if(!checkFasta(originalBuild)|!checkFasta(newBuild)) stop("Error: originalBuild and/or newBuild do not appear to be valid fasta files")
  if(is.null(tmpdir))
  {
    tmpdir = tempdir()
  }
  if(!dir.exists(tmpdir)) tryCatch(dir.create(tmpdir),error = function(e){print("Could not create temp directory at supplied path. see below for details:");stop(e)})

  indexFasta(originalBuild)
  indexFasta(newBuild)
  originalFa = FaFile(file = originalBuild)
  seqinfo(originalFa)@seqnames
  unique(chrom)
  snpsSeq = tryCatch(getSeq(x=originalFa,GRanges(seqnames=chrom,ranges = IRanges(start=start-lengthSides,end=end+lengthSides))),error=function(e){print("failed to find locations specified in fasta file. see below for details:");stop(e)})
  if(is.null(varnames))
  {
    varnames = paste(chrom,start,end,sep="_")
  }
  names(snpsSeq) = varnames
  writeXStringSet(x = snpsSeq,format="fastq",filepath=sprintf("%s/liftOverInput.fastq",tmpdir))
  write.table(data.frame(FileName=sprintf("%s/liftOverInput.fastq",tmpdir),SampleName=sprintf("%s/liftOverInput",tmpdir)),file=sprintf("%s/liftOverInput",tmpdir),quote=F,sep="\t",col.names=T,row.names=F)
  aligned = qAlign(sprintf("%s/liftOverInput",tmpdir),newBuild ,cacheDir = tmpdir,maxHits = 1)

  algn = scanBam(aligned@alignments$FileName)
  newInOld = match(algn[[1]]$qname,varnames)
  newRegs = data.frame(name=algn[[1]]$qname,
                       chrom_orig = chrom[newInOld],
                       start_orig = start[newInOld],
                       end_orig = end[newInOld],
                       chrom = NA,
                       start = NA,
                       end = NA,
                       stringsAsFactors = F)

  choos = which(algn[[1]]$cigar==sprintf("%sM",2*lengthSides+1)&!is.na(algn[[1]]$rname))
  newRegs[choos,"chrom"] = as.character(algn[[1]]$rname[choos])
  newRegs[choos,"start"] = algn[[1]]$pos[choos]+lengthSides
  newRegs[choos,"end"] = algn[[1]]$pos[choos]+lengthSides + end[choos]-start[choos]
  newRegs
}
