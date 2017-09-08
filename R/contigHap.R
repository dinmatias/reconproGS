#' Assembles the contig haplotype from a phased vcf
#'
#' Merge SNPs from the same contig to generate haplotype data; the vcf
#' file must be phased.
#'
#' @param filePath full path to the vcf file; this calls the vcfLoad
#' @param dataIn data.table obtained using vcfLoad function
#' @return A data.table with the contig name, sample id and haplotype.
#' @examples
#' # To follow


contigHap <- function(filePath = NULL, dataIn){

  if(!is.null(filePath)){
    dataIn <- vcfLoad(filePath = filePath)
  }

  if("|" %in% unlist(strsplit(dataIn$GT[1], ""))){

    sepPhase <- "|"

  } else if("/" %in% unlist(strsplit(dataIn$GT[1], ""))){

    print("The VCF is not PHASED")
    stop()

  }else{

    print("Cannot identify if VCF is PHASED or NOT")
    stop()

  }

  subData <- recodeGT(dataIn)

  x <- names(subData)[grep(pattern = "Allele", x = names(subData))]

  commandGen <- parse(text = paste(".(",
                                   paste(paste('Hap', 1:length(x),
                                               ' = paste(c(Allele' ,
                                               1:length(x),
                                               '), collapse = "")',
                                               sep = ""),
                                         collapse = ","),
                                   ")", sep = ""))

  newData <- subData[ , eval(commandGen), by = .(`#CHROM`,SampleID)]

  return(newData)

}

