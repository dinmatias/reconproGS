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

getAlleleCount <- function(tempRes){

  tempRes[ , identifier := 1]

  tempSub <- data.table::dcast(tempRes,  `#CHROM` + POS + V1 ~ SampleID ,
                               fun.aggregate = length, value.var = "identifier")

  namesSamp <-unique(tempRes$SampleID)

  partialRes <- data.table::melt(tempSub,
                                 id.vars = colnames(tempSub)[
                                   !(colnames(tempSub) %in% namesSamp)],
                                 measure.vars = namesSamp)

  makeNA <- partialRes[is.na(V1) & (value > 0), ]

  partialRes[`#CHROM` %in% makeNA$`#CHROM` &
              POS %in% makeNA$POS &
              variable %in% makeNA$variable, value := NA]

  partialRes <- partialRes[!is.na(V1), ]

  names(partialRes)[3:5] <- c("Allele", "SampleID", "Count")

  return(partialRes)

}
