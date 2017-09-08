#' Recode the genotype with the allele
#'
#' Converts the genotype written by vcf using numbers to allele. This is
#' internally called by other functions that are involved in getting allele
#' counts. Moreover, this functions simplifies the data by removing other
#' information of the vcf.
#'
#' @param dataIn data.table obtained using vcfLoad function
#' @return A data.table with contig info (#CHROM), SNP position (POS), sample ID
#' and the alleles of the sample for that particular SNP
#' @examples
#' # To follow

recodeGT <- function(dataIn){

  sepPhase <- unlist(strsplit(dataIn$GT[1], ""))[[2]]

  ploidy <- length(unlist(strsplit(dataIn$GT[1], sepPhase, fixed = TRUE)))

  hapNames <- paste("Hap", 1:ploidy, sep = "")

  dataIn[ , (hapNames) :=
            data.table::tstrsplit(GT, sepPhase, fixed = TRUE)][
              , (hapNames) := lapply(.(Hap1, Hap2), as.numeric)
              ]

  dataIn[!(is.na(GT)), paste("ALT", 1:ncol(dataIn[ , data.table::tstrsplit(ALT, ",", fixed = TRUE)]),
                  sep = "") := data.table::tstrsplit(ALT, ",", fixed = TRUE)]


  subData <- dataIn[ , c("#CHROM", "POS", "SampleID", hapNames, "REF",
                         names(dataIn)[grep(x = names(dataIn),
                                            pattern = "ALT[0-9].*$")]),
                     with = FALSE]

  rm(dataIn)

  alleleNames <- paste("Allele", 1:ploidy, sep = "")

  tempMark <- length(c("#CHROM", "POS", "SampleID", hapNames, "REF"))

  subData[ , (alleleNames)] <- lapply(subData[ , (hapNames), with = FALSE],
                                                function(x) as.matrix(subData[ , tempMark:ncol(subData)]
                                                )[cbind(seq_len(nrow(subData)), x + 1)])

  x <- c("#CHROM", "POS", "SampleID", alleleNames)

  subData <- subData[ , x, with = FALSE]

  return(subData)

}
