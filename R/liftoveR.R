#' @import Rbowtie
#' @import Rsamtools
#' @import Biostrings
#' @import GenomicRanges
#' @import stringr
#' @import IRanges
#' @importFrom utils read.table
#' @importFrom utils write.table

checkBowtiePrefix = function(pref)
{
  if(file.exists(sprintf("%s.1.bt2")))
    return(T)
  return(F)
}

indexFasta = function(fa)
{
  sprintf("Building fasta index for %s",fa)
  if(!file.exists(sprintf("%s.fai",fa)))
    indexFa(fa)
}

subtractSides = function(pos,lengthSides)
{
  lengthSides = rep(lengthSides,length(pos))
  lengthSides[(pos-lengthSides)<1] = pos[(pos-lengthSides)<1] - 1
  lengthSides
}

#' Register paths to blast and samtools programs
#' @param blastn path to the blastn executable. Currently, only version 2.2.31+ supported, since it's the only one with appropriate SAM output
#' @param makeblastdb path to the makeblastdb executable.
#' @param samtools path to the samtools executable
#'@export
registerBlast = function(blastn,makeblastdb,samtools)
{
  if(!file.exists(blastn) | !file.exists(makeblastdb) | !file.exists(samtools))
    stop("Error: invalid path given to blastn, makeblastdb or samtools")
  versionStringBlast = system(sprintf("%s -version",blastn),intern=T)[1]
  versionStringBlastDB = system(sprintf("%s -version",makeblastdb),intern=T)[1]
  versionStringSamtools = system(sprintf("%s --version",samtools),intern=T)[1]
  if(!grepl("blastn",versionStringBlast))
    stop("Error: incorrect file referenced as blastn")
  if(!grepl("makeblastdb",versionStringBlastDB))
    stop("Error: incorrect file referenced as makeblastdb")
  if(!grepl("samtools",versionStringSamtools))
    stop("Error: incorrect file referenced as samtools")
  assign("blastn",blastn,envir=.blastEnv)
  assign("makeblastdb",makeblastdb,envir=.blastEnv)
  assign("samtools",samtools,envir=.blastEnv)
  print("blast programs registered")
}

findDisplacementByCigar = function(cigar,lengthSide,side)
{
  lengthSide = as.numeric(lengthSide)
  if(is.na(cigar))
    return(NA)
  if(!side %in% c("start","end"))
    stop("side is either start or end")
  numbers = stringr::str_split(cigar,"[aA-zZ]+",simplify = T)[1,]
  numbers = as.numeric(numbers[-length(numbers)])
  types = stringr::str_split(cigar,"[0-9]+",simplify = T)[1,-1]
  if(side == "end")
  {
    displace = sum(numbers[types!="D"]) - lengthSide - 1
    pos = which(types!="I")[which(cumsum(numbers[types!="I"]) > displace)[1]]

  } else
  {
    displace = lengthSide
    pos = which(types!="I")[which(cumsum(numbers[types!="I"]) > displace)[1]]
  }
  if(types[1]=="S"|types[1]=="H")
    displace = displace - numbers[1]
  if(any(types[1:pos] == "D"))
    displace = displace + sum(numbers[1:pos][types[1:pos]=="D"])
  if(any(types[1:pos] == "I"))
    displace = displace - sum(numbers[1:pos][types[1:pos]=="I"])
  return(displace)
}

addSides = function(chrom,pos,lengthSides,originalBuild)
{
  seqTable = data.frame(scanFaIndex(originalBuild))
  chromLengths = seqTable$width[match(chrom,seqTable$seqnames)]

  lengthSidesVector = rep(lengthSides,length(pos))
  toChange = (pos+lengthSidesVector)>chromLengths
  if(any(toChange))
    lengthSidesVector[toChange] = chromLengths[toChange] - pos[toChange]
  lengthSidesVector
}

checkFasta = function(fa)
{
  firstLine = read.table(fa,nrows=1,stringsAsFactors = F)
  if(!substr(firstLine,1,1)[1] == ">")
    return(F)
  T
}

getSequences =  function(originalFa,chrom,start,end)
{
  tryCatch(getSeq(x=originalFa,GRanges(seqnames=chrom,
                                       ranges = IRanges(start=start,
                                                        end))),
           error=function(e){print("failed to find locations specified in fasta file. see below for details:");stop(e)})
}

bowtieAlign = function(reads,ref,cacheDir)
{
  RBowtieRefDir = sprintf("%s/%s.RBowtie",cacheDir,basename(ref))
  if(!file.exists(paste(RBowtieRefDir,"/index.1.ebwt",sep="")))
  {
    print(sprintf("Index file for RBowtie not found, building in directory %s/%s.RBowtie",cacheDir,basename(ref)))
    Rbowtie::bowtie_build(references = ref,outdir = RBowtieRefDir,force = T)
  }
  outFile = paste(reads,".sam",sep="")
  print("Performing Bowtie alignment")
  Rbowtie::bowtie(sequences = reads,index = paste(RBowtieRefDir,"/index",sep=""),type = "single",force = T,outfile = outFile,S=T)
  bamFile = Rsamtools::asBam(outFile,outFile)
  sortedBamFile = Rsamtools::sortBam(bamFile,paste(bamFile,"sorted",sep=""))
  return(sortedBamFile)
}

doAlignment = function(seqs,ref,tmpdir,aligner,memlimit,maxMismatches,blastArgs = "")
{
  if(aligner == "Rbowtie")
  {
    aligned = bowtieAlign(seqs,ref ,cacheDir = tmpdir)
    return(aligned)
  } else if(aligner =="Rsubread")
  {
    if(!file.exists(sprintf("%s.reads",ref)))
    {
      print(sprintf("Index file for Rsubread not found, building index in %s.*",ref))
      Rsubread::buildindex(basename = ref,reference = ref,memory = memlimit)
    }
    Rsubread::align(index =ref,readfile1 = sprintf("%s",seqs),type = 1,phredOffset = 64,maxMismatches = maxMismatches,complexIndels=T,indels=200,minFragLength = 25)
    return(sprintf("%s.subread.BAM",seqs))
  } else if(aligner == "blastn")
  {
    blastn = get("blastn",envir = .blastEnv)
    makeblastdb = get("makeblastdb",envir = .blastEnv)
    samtools = get("samtools",envir = .blastEnv)
    if(is.null(blastn)|is.null(makeblastdb)|is.null(samtools))
    {
      stop("blastn, makeblastdb or samtools not registered")
    }
    blastdb = file.path(tmpdir,"blastdb")
    bampath = file.path(sprintf("%s.blast.bam",seqs))
    sampath = file.path(sprintf("%s.blast.sam",seqs))
    print("building blast database")
    system(sprintf("%s -in %s -out %s -dbtype nucl -parse_seqids",makeblastdb,ref,blastdb))
    print("performing alignment (may take a long time, this IS blast)....")
    system(sprintf("%s -db %s -query %s %s -out %s",blastn,blastdb,seqs,blastArgs,sampath),ignore.stderr = T,ignore.stdout = T)
    system(sprintf("%s view -bS %s -o %s",samtools,sampath,bampath),ignore.stderr = T,ignore.stdout = T)
    return(bampath)
    }
}

#' Lift over coordinates of variants from one genome build to another, by realignment
#' @param chrom For each variant, the chromosome on which it resides
#' @param start,end Physical coordinates of variants (in 1-based indexing. see \url{http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms})
#' @param originalBuild Path to the fasta sequence of the build from which the variants originate
#' @param newBuild Path to the fasta sequence of the build to which the variants will be lifted over
#' @param varnames Variant names. If missing, variants will be given names based as chromosome_position in original build
#' @param lengthSides Number of bases which will be extracted from each side of the variant
#' @param maxMismatches Corresponds to the number of mismatches. See the vignette for more information on how this is specified for each aligner. {Not implemented for blastn}
#' @param tmpdir Working directory for sequence alignment etc. Defaults to a temporary directory created automatically by the current R session.
#' @param aligner use Rbowtie, Rsubread or blastn for alignment? Recommended that be set to blastn, but could be slow on large genomes, and unneccessary for aligning simple polymorphisms.
#' @param memlimit Memory (in MB) limit for building the Rsubread index (relevant only if Rsubread is TRUE). Defaults to 8000MB, which should allow mammalian genomes to be built with high efficiency for fast alignment (see Rsubread documentation for more detail)
#' @param maxAllowedLength the maximum length of read that will be aligned (including the sides). regions longer than that will be split and rejoined after alignment
#' @param blastArgs advanced: customize the arguments to be passed to blast. by default, blast will be run with with "-max_hsps 1 -max_target_seqs 1 -outfmt 15 -parse_deflines". These are required for the parsing later to work, so please add those in addition to additional argument desired. (user responsibility - no validation by the package)
#' @return A data frame containing the variants with their positions in the new build, based on the following scheme:
#' \describe{
#' \item{name}{Name of the variant}
#' \item{chrom_orig, start_orig, end_orig}{The original position supplied to function.}
#' \item{chrom, start, end}{Coordinates of alignment in the new build}
#' }
#' @export
liftover = function(chrom,start,end = NULL,originalBuild,newBuild,varnames=NULL,lengthSides=25,maxMismatches = 0,tmpdir = NULL,aligner = "Rsubread",memlimit = 8000, maxAllowedLength = 500,blastArgs= "-max_hsps 1 -parse_deflines -outfmt 15 -max_target_seqs 1")
{
  if(aligner == "Rsubread")
  {
    if(!requireNamespace("Rsubread", quietly = TRUE))
    {
      stop("Rsubread is not installed. Please install from Bioconductor. Note: Rsubread only supports Linux and OSX")
    }
  }
  if(!is.numeric(lengthSides)|length(lengthSides)!=1|lengthSides<0|round(lengthSides)!=lengthSides) stop("Error: lengthSides must be a single positive round number")
  if(is.null(end)) end = start
  if(!is.numeric(start)|!is.numeric(end)) stop("Error: Start and/or End positions supplied are not numeric vectors!")
  rangeCheck = range(length(start),length(end),length(chrom))
  if(!rangeCheck[1]==rangeCheck[2]) stop("Error: chrom, start and end (if supplied) MUST have the same length")
  if(!all(end>=start)) stop("Error: all End positions must be at >= Start positions")
  if(!file.exists(originalBuild)|!file.exists(newBuild)) stop("Error: originalBuild and/or newBuild are not paths to existing files")
  if(!checkFasta(originalBuild)|!checkFasta(newBuild)) stop("Error: originalBuild and/or newBuild do not appear to be valid fasta files")
  if(is.null(tmpdir))
  {
    tmpdir = tempdir()
  }
  if(is.null(varnames))
  {
    varnames = paste(1:length(start),chrom,start,end,sep="_")
  }
  if(length(varnames)!=length(start)) stop("Error: if supplied varnames MUST be be same length as chrom and start")
  if(!dir.exists(tmpdir)) tryCatch(dir.create(tmpdir),error = function(e){print("Could not create temp directory at supplied path. see below for details:");stop(e)})
  if(aligner == "blastn")
  {
    if(is.null(.blastEnv$blastn)|is.null(.blastEnv$makeblastdb)|is.null(.blastEnv$samtools))
    {
      stop("blastn, makeblastdb or samtools not registered")
    }
    readsPath = file.path(tmpdir,"liftOverInput.fasta")
    readsFormat = "fasta"
  } else
  {
    readsPath = file.path(tmpdir,"liftOverInput.fastq")
    readsFormat = "fastq"
  }

  newRegsCollapsed = data.frame(name=varnames,chrom_orig=chrom,start_orig=start,end_orig=end,chrom=NA,start=NA,end=NA,stringsAsFactors = F)


  chrom = as.character(chrom)
  indexFasta(originalBuild)
  indexFasta(newBuild)
  originalFa = FaFile(file = originalBuild)
  seqLengths = end - start + 1
  nsplits = floor(seqLengths/(maxAllowedLength))
  print(sprintf("Identified %s regions that will need to be split due to size",sum(nsplits>0)))
  cat(sprintf("Splitting Reads..."))
  whichSplit = which(nsplits>0)
  if(length(whichSplit)>0)
  {
      if(lengthSides<25&aligner=="Rsubread")
        stop("Selected Rsubread aligner and regions need to be split. lengthSides MUST be at least 25 to ensure split reads are long enough for Rsubread to align!")
      k = 0
      splitRegs = NULL
      splitStart = NULL
      splitEnd = NULL
      splitChrom = NULL
      splitNames = NULL
      for(toSplit in whichSplit)
      {
        k = k+1
        if(k %% floor(length(whichSplit)/10)==0)
          cat(sprintf("..%s%%",ceiling(k/length(whichSplit)*100)))
        for(split in 1:(nsplits[toSplit]+1))
        {
          curStart = start[toSplit] + (split-1)*maxAllowedLength
          splitStart = c(splitStart,curStart)
          splitEnd = c(splitEnd,min(end[toSplit],curStart + maxAllowedLength - 1))
          splitChrom = c(splitChrom,chrom[toSplit])
          splitNames = c(splitNames,paste(varnames[toSplit],split,sep="splitNumber"))
        }
      }
      start = c(start[-whichSplit],splitStart)
      end = c(end[-whichSplit],splitEnd)
      chrom = c(chrom[-whichSplit],splitChrom)
      varnames = c(varnames[-whichSplit],splitNames)
  }

  leftSide = subtractSides(start,lengthSides)
  rightSide = addSides(chrom,end,lengthSides,originalBuild)
  snpsSeq = getSequences(originalFa,
                         chrom,
                         start - leftSide,
                         end + rightSide)

  names(snpsSeq) = varnames
  writeXStringSet(x = snpsSeq,format=readsFormat,filepath=readsPath)

  aligned = doAlignment(seqs = readsPath,ref = newBuild ,tmpdir = tmpdir,aligner = aligner,memlimit = memlimit,maxMismatches = maxMismatches,blastArgs)
  param <- ScanBamParam(what=c("rname","pos","qname","cigar","strand"), tag=c("MD", "NM"))
  algn = scanBam(aligned,param = param)
  if(aligner == "blastn")
  {
    algn[[1]]$qname = sub("lcl\\|","",algn[[1]]$qname)
    algn[[1]]$rname = sub("lcl\\|","",algn[[1]]$rname)

  }
  newInOld = match(algn[[1]]$rname,varnames)
  newRegs = data.frame(name=algn[[1]]$qname,
                       chrom_orig = chrom[newInOld],
                       start_orig = start[newInOld],
                       end_orig = end[newInOld],
                       chrom = NA,
                       start = NA,
                       end = NA,
                       stringsAsFactors = F)
  choos = which(!is.na(algn[[1]]$rname))
  if(aligner == "Rbowtie")
  {
    choos = intersect(which(algn[[1]]$tag$NM<=maxMismatches),choos)
  }
  print("Parsing new start positions based on alignment")
  newRegs[choos,"start"] = algn[[1]]$pos[choos]+ apply(cbind(algn[[1]]$cigar[choos],leftSide[choos]),1,function(x)
                                                      findDisplacementByCigar(x[1],x[2],side = "start"))
  print("Parsing new end positions based on alignment")
  newRegs[choos,"end"] = algn[[1]]$pos[choos]+apply(cbind(algn[[1]]$cigar[choos],rightSide[choos]),1,function(x)
                                                        findDisplacementByCigar(x[1],x[2],side = "end"))
  newRegs[choos,"chrom"] = as.character(algn[[1]]$rname[choos])
  groups = gsub("splitNumber+.*","",newRegs$name)

  finalStarts = tapply(newRegs$start,groups,min)
  finalEnd = tapply(newRegs$end,groups,max)
  finalChrom = tapply(newRegs$chrom,groups,function(x)ifelse(length(unique(x))==1,x[1],NA))
  newRegsCollapsed$start[match(names(finalStarts),newRegsCollapsed$name)] = finalStarts
  newRegsCollapsed$end[match(names(finalEnd),newRegsCollapsed$name)] = finalEnd
  newRegsCollapsed$chrom[match(names(finalChrom),newRegsCollapsed$name)] = finalChrom
  anyIsNA = apply(newRegsCollapsed[,5:7],1,function(x)any(is.na(x)))
  newRegsCollapsed[anyIsNA,5:7] = NA
  newRegsCollapsed
}
