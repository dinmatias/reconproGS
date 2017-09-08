#' Fast conversion of vcf to genind
#'
#' This function converts a raw vcf file to a genind object with the
#' assumptions that (1) the ploidy is 2 (but will work with one actually,
#' but not with ploidy > 2) and (2) marker is codominant (duhhhhhh, it
#' should be!!!) this basically creates an allele frequency table from
#' the vcf which is the used as the tab in a genind object. Then, the
#' output is wrapped using the genind function from adegenet together
#' with other auxillary information.
#'
#' @param filePath full path to the vcf file; this calls the vcfLoad
#' @param dataIn data.table obtained using vcfLoad function
#' @param pop vector of the the popID of each sample
#' @param ... other parameters needed by genind type object, see genind or df2genind
#' @return A genind class object to be used by adegenet
#' @examples
#' # reads the example.vcf and makes a genind object
#' vcfGenind <- vcf2genind("example.vcf")
#' populationData <- read.table("populationMetadata")
#' findfixedSNP(dataIn = vcfInput, popinfo = populationData)

vcf2genind <- function(filePath = NULL, inputData = NULL, pop = NULL, prevcall = NULL,
                       ploidy = NULL, type = "codom",
                       strata = NULL, hierarchy = NULL){

  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop("Pkg needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # reads the vcf file
  dataIn <- vcfLoad(filePath = filePath)

  if(is.null(ploidy)){

    sepPhase <- unlist(strsplit(dataIn[!is.na(GT), GT][1], ""))[[2]]

    ploidy <- length(unlist(strsplit(dataIn[!is.na(GT), GT][1], sepPhase, fixed = TRUE)))
  }

  tempRes <- recodeGT(dataIn)

  x <- length(names(tempRes)) - 3

  commandGen <- parse(text = paste('c(', paste(paste("Allele", 1:x, sep = ""),
                                               collapse = ","),
                                   ')'))

  tempRes <- tempRes[ , eval(commandGen), by = .(SampleID, `#CHROM`, POS)]

  partialRes <- getAlleleCount(tempRes)

  # removes tempRes and collect trash
  rm(tempRes)

  partialRes[ , markerID := paste(`#CHROM`, POS, sep = "__")][
    , `#CHROM` := NULL][ , POS := NULL]

  partialRes[ , AlleleID := paste(markerID, Allele, sep =".")][
    , markerID := NULL]

  partialRes <- dcast(partialRes,  SampleID ~ AlleleID, value.var = "Count")

  outFile <- as.data.frame(partialRes[ , -1])

  row.names(outFile) <- as.character(partialRes[ , SampleID])

  rm(partialRes)

  # make a genind object
  out <- adegenet::genind(tab = outFile, pop = pop, prevcall = prevcall,
                          ploidy = ploidy, type = type, strata = strata, hierarchy = hierarchy)

  # return a genind object
  return(out)

}

