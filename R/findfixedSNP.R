#' Calculates allele frequency in each population.
#'
#' Determines the allele frequencies of each SNPs in specifed populations
#'
#' @param filePath full path to the vcf file; this calls the vcfLoad
#' @param dataIn data.table obtained using vcfLoad function
#' @param popinfo dataframe with  population identifier (column1) and sample identifier (column2)
#' @return A data.table with population.allele as columns and SNPs as rows with frequencies as values
#' @examples
#' # reads example.vf, reads table containing sample id and population
#' # id of each sample then computes allele frequency for each snp in
#' # each population
#' vcfInput <- vcfLoad("example.vcf")
#' populationData <- read.table("populationMetadata")
#' findfixedSNP(dataIn = vcfInput, popinfo = populationData)

findfixedSNP <- function(filePath = NULL, dataIn, popinfo){

  if(!is.null(filePath)){
    dataIn <- vcfLoad(filePath = filePath)
  }

  # make sure data.table format, and names are correct
  popinfo <- data.table(popinfo)
  names(popinfo) <- c("Pop", "Samples")
  popinfo$Pop <- as.character(popinfo$Pop)
  popinfo$Samples <- as.character(popinfo$Samples)

  # use sampleID as key
  setkey(dataIn, SampleID)
  setkey(popinfo, Samples)

  # include the population info to the data table
  dataIn <- dataIn[popinfo]

  # generate a locus SNP identifier; this one is easier
  # than using two columns (will results to 3D matrix rather than 4D)
  dataIn[ , locusID := paste(`#CHROM`, POS, sep = "____") ]
  setkey(dataIn, locusID)

  # use the genotype information as counter for alternative allele
  # since 0 in the genotype indicate reference allele and 1 indicate alternative
  # this parse the genotype and convert then to numeric
  dataIn[ , c("freqA1", "freqA2") := tstrsplit(dataIn$GT, "/") ][
    , c("freqA1", "freqA2") := lapply(.(freqA1, freqA2), as.numeric)]

  # split the ALT allele
  dataIn[ , paste("ALT", 1:ncol(dataIn[ , tstrsplit(ALT, ",")]), sep = "") := tstrsplit(ALT, ",") ]

  # just get the necessary columns and ditch the other crap
  # retain locusID, sampleID and the alleles for that sample
  subData <- dataIn[ , c("locusID", "SampleID", "Pop", "freqA1",
                         "freqA2", "REF",
                         names(dataIn)[grep(x = names(dataIn),
                                            pattern = "ALT[0-9].*$")]),
                     with = FALSE]

  # Recode the allele freq with the "base" given in the REF/ALT columns
  # this is hard-coded for up to two alternative alleles only
  # need to change this one if there are more OR write a code that identifies the number
  # of alleles and change accordingly
  subData[ , c("Allele1", "Allele2")] <- lapply(subData[ , c("freqA1", "freqA2"), with = FALSE],
                                                function(x) as.matrix(subData[ , 6:ncol(subData)]
                                                )[cbind(seq_len(nrow(subData)), x + 1)])


  # concatenates all alleles in a population per locus
  # output is a vector of all alleles per population per SNP
  tempRes <- subData[ , c(Allele1, Allele2), by = .(Pop, locusID)]

  # makes a 3D matrix making getting the frequency of allele
  # per in each locus per population
  # results[allele, population, locusID]
  tempRes <- table(tempRes$V1, tempRes$Pop, tempRes$locusID)

  # convert the 3D matrix into a list of 2D matrix
  partialRes <- lapply(seq(dim(tempRes)[3]), function(x) tempRes[,,x])

  # rename the lsit with locusName
  names(partialRes) <- dimnames(tempRes)[[3]]

  # computes the allelic frequency in each population
  # then converts each matrix into a data.table
  freqData <- lapply(partialRes, function(x) as.data.table(apply(x, 2, function(d) d/sum(d)),
                                                           keep.rownames = TRUE) )

  # makes one master data.table
  # with the locus name in .id column and alleles in rn colum
  freqData <- rbindlist(freqData, idcol = TRUE)

  # removes alleles that are not present
  outFile <- freqData[-which(apply(freqData[ , 3:ncol(freqData)], 1, sum) == 0), ]

  # determines the private allele
  # to be a private allele, it must be present in at least one population
  # AND it should be zero in atleast one population
  privAlleles <- apply(outFile[ , 3:ncol(outFile)], 1,  function(x) ifelse( ((sum(x == 0) >= 1) & (sum(x > 0) >= 1)), 1, 0))

  # determines the fixed alleles
  # to be a fixed allele, it must have a frequency of 1 in least one population
  # AND it should be zero in atleast one population
  fixAlleles <- apply(outFile[ , 3:ncol(outFile)], 1,  function(x) ifelse( ((sum(x == 0) >= 1) & (sum(x == 1) >= 1)), 1, 0))

  outFile <- cbind(outFile, privAlleles, fixAlleles)


  # separates locus name and pos to two different columns
  outFile[, c("#CHROM", "POS") := tstrsplit(.id, "____", fixed=TRUE)][ , .id := NULL]

  # reorder columns
  setcolorder(outFile, colnames(outFile)[c(
    (ncol(outFile) - 1), ncol(outFile), 1:(ncol(outFile)-2))])

  # return the three data tables
  return(outFile)

}
