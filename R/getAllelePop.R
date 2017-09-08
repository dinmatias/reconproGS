#' Determines the number of copies of each allele of a locus in every sample
#'
#' Counts the number of each allele of a locus in each individual. If an individual
#' or sample is not genotyped for a certain locus, the count will be marked as NA.
#'
#' @param tempRes data.table summarizing the alleles present in all samples.
#' This should be an output of recodeGT function.
#' @return A data.table with contig info (#CHROM), SNP position (POS), Allele,
#' sample ID, and the count of the allele of a particular SNP in an individual
#'
#' @examples
#' # To follow

getAllelePop <- function(dataIn, popInfo, ploidy = 2){


  dataIn <- recodeGT(dataIn)

  x <- length(names(dataIn)) - 3

  commandGen <- parse(text = paste('c(', paste(paste("Allele", 1:x, sep = ""),
                                               collapse = ","),
                                   ')'))

  dataIn <- dataIn[ , eval(commandGen), by = .(SampleID, `#CHROM`, POS)]

  dataIn <- getAlleleCount(dataIn)

  setkey(dataIn, SampleID)

  dataIn <- dataIn[popInfo]

  tempOut <- data.table::dcast(dataIn,  `#CHROM` + POS + Allele ~ Pop ,
                               fun.aggregate = sum, value.var = "Count")


  namesPop <-unique(popInfo$Pop)

  tempOut <- data.table::melt(tempOut,
                              id.vars = colnames(tempOut)[
                                !(colnames(tempOut) %in% namesPop)],
                              measure.vars = namesPop)

  return(tempOut)

}
