#' Determines the frequency of each allele in each population
#'
#' Takes input vcf or data.table loaded using vcfLoad function and determines
#' the allelic frequency of every allele in each SNP in the provided population
#' information
#'
#' @param filePath full path to the vcf file
#' @param dataIn data.table obtained using vcfLoad function
#' @param popInfo dataframe with  population identifier (column1) and sample
#' identifier (column2)
#' @return A data.table with contig info (#CHROM), SNP position (POS), Allele,
#' Population ID (Pop), and the count of the allele (Count), estimate of allele
#' freuquency in a population (Frequency), and variance of the allele frequency
#' (Variance).
#'
#' @examples
#' # To follow

alleleFrequency <- function(filePath = NULL, dataIn = NULL, popInfo = NULL, ploidy = 2){

  if(is.null(filePath) & is.null(dataIn)){
    stop("Provide a file path to a vcf file or a data.table loaded by vcfLoad")
  }

  if(is.null(popInfo)){
    stop("Population information must be given")
  }

  if(!is.null(filePath)){
    dataIn <- vcfLoad(filePath)
  }

  popInfo <- data.table(popInfo)

  names(popInfo) <- c("Pop", "Samples")

  popInfo[ , .(Pop = as.character(Pop), Samples = as.character(Samples))]

  setkey(popInfo, Samples)


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


  tempOut[ , Frequency := value/sum(value, na.rm = TRUE),
           by = .(`#CHROM`, POS, variable)]

  names(tempOut)[3:5] <- c("Allele", "Pop", "Count")

  numGen <- dataIn[ , .(numGen = sum(Count, na.rm = TRUE)/ploidy),
                     by = .(`#CHROM`, POS, Pop)]

  setkey(numGen, `#CHROM`, POS, Pop)

  numHomo <- dataIn[ , .(homo = length(which(Count == 2))),
                     by = .(`#CHROM`, POS, Allele, Pop)]

  setkey(numHomo, `#CHROM`, POS, Pop)

  homozygous <- numHomo[numGen]

  setkey(homozygous, `#CHROM`, POS, Allele, Pop)

  setkey(tempOut, `#CHROM`, POS, Allele, Pop)

  tempOut <- tempOut[homozygous]

  # Genetic Data Analysis II Weir 1996 p 39
  tempOut[ , Variance := (Frequency + (homo/numGen) - (2 * Frequency^2))/(ploidy * numGen)]

  tempOut[ , c("homo", "numGen") := NULL]

  setkey(tempOut, Pop, `#CHROM`, POS, Allele)

  return(tempOut)

}
