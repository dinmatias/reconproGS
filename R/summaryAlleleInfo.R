summaryAlleleInfo <- function(dataIn){

  # generate a locus SNP identifier; this one is easier
  # than using two columns (will results to 3D matrix rather than 4D)
  dataIn[ , locusID := paste(`#CHROM`, POS, sep = "____") ]
  setkey(dataIn, locusID)

  # use the genotype information as counter for alternative allele
  # since 0 in the genotype indicate reference allele and 1 indicate alternative
  # this parse the genotype and convert then to numeric
  dataIn[ , c("freqA1", "freqA2") := tstrsplit(dataIn$GT, "/") ][
    , c("freqA1", "freqA2") := lapply(.(freqA1, freqA2), as.numeric)]

  # works with 3 possible alleles at the moment
  dataIn[!(is.na(GT)) , paste("ALT",
                              1:ncol(dataIn[ , tstrsplit(ALT, ",")]),
                              sep = "") := tstrsplit(ALT, ",") ]

  # just get the necessary columns and ditch the other crap
  # retain locusID, sampleID and the alleles for that sample
  subData <- dataIn[ , c("locusID", "SampleID", "freqA1",
                         "freqA2", "REF",
                         names(dataIn)[grep(x = names(dataIn),
                                            pattern = "ALT[0-9].*$")]),
                     with = FALSE]

  # recode the genotype with the allele
  subData[ , c("Allele1", "Allele2")] <- lapply(subData[ , c("freqA1", "freqA2"), with = FALSE],
                                                function(x) as.matrix(subData[ , 5:ncol(subData)]
                                                )[cbind(seq_len(nrow(subData)), x + 1)])

  # just get the necessary columns and ditch the other crap
  # retain locusID, sampleID and the alleles for that sample
  subData <- subData[ , .(locusID, SampleID, Allele1, Allele2) ]

  # concatenates the alleles per locus per individual
  tempRes <- subData[ , c(Allele1, Allele2), by = .(SampleID, locusID)]

  d <- x[ , c(Allele1, Allele2), by = .(SampleID, `#CHROM`, POS)]

  return(tempRes)
}
