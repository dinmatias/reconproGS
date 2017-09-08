#' Loads a vcf and parse by SNP information per sample
#'
#' @param filePath full path to the vcf file
#' @param SNPstart number specifying which SNP to start to read
#' @param SNPend number specifying which SNP to read last
#' @return A data.table of parsed SNP in long format
#' @examples
#' # reads all SNPs from the example.vcf
#' vcfLoad(filePath = "example.vcf")
#' # reads the 500th up to 3000th SNP
#' vcfLoad(filePath = "example.vcf", SNPstart = 500, SNPend = 3000)

vcfLoad <- function(filePath = NULL, SNPstart = NA, SNPend = NA){

  # control flow for filePath
  if(is.null(filePath)){
    stop("filePath must be specified")
  }

  if(!is.na(SNPstart) & !is.na(SNPend)){
    if(!is.numeric(SNPend) | SNPend%%1 > 0){
      stop("Specify proper last SNP NUMBER to read in SNPstart" )
    }
    if(!is.numeric(SNPstart) | SNPstart%%1 > 0){
      stop("Specify proper starting SNP NUMBER to read in SNPstart" )
    }
    # control for SNP subset
    if(SNPstart > SNPend){
      stop("SNPend must be greater than SNPstart")
    }

  }

  fileData <- data.table::fread(filePath, header = T, stringsAsFactors = F,
                                na.strings = c("","NA","#N/A"), fill = T,
                                blank.lines.skip = FALSE,
                                skip= "#CHROM")

  # subset the file data IF SNPstart and SNPend were specified
  if(!is.na(SNPstart) & !is.na(SNPend)){
    # subset the data.table from rows SNPstart to SNPend
    fileData <- fileData[SNPstart:SNPend, ]
  }

  # get the names of additional information placed in the INFO column of vcf;
  # each variable is separated by ";"; "=.*$ is used to remove assigned value
  namesInfo <- gsub(x = strsplit(x = fileData$INFO[1], split = ";")[[1]],
                    pattern = "=.*$", replacement = "")

  # determine the information associated with the genotype for each invidual
  # SNP will be used as column names when genotype information is parsed
  # each variable is separated by ":"
  namesForm <- strsplit(x = fileData$FORMAT[1], split = ":")[[1]]

  # determine the sample names (column names of genotype fields) which are
  # the columns after the genotype format field (FORMAT column of vcf)
  # possible to change this and make use of the fact that there are 9 colums
  # prior to genotype fields if the vcf has genotype information
  namesSamp <- colnames(fileData)[10:ncol(fileData)]

  # each sample has it's own column; here, it will be melted to 2 columns, one the sample
  # identifier (which is the sample name) and the sample info (the values in the species column)
  # melt sample columns and make the sample name a variable and the sample info as value
  reformDat <- data.table::melt(fileData,
                                id.vars = colnames(fileData)[!(colnames(fileData) %in% namesSamp)],
                                measure.vars = namesSamp)

  # remove dataframe after converting to reformDat
  rm(fileData)

  # rename the melted variable and convert to strings
  names(reformDat)[names(reformDat) == "variable"] <- "SampleID"
  reformDat$SampleID <- as.character(reformDat$SampleID)

  # split the info about each SNP to different columns
  # the content of this column is specified in the namesForm
  if(length(namesForm) == 1){
    names(reformDat)[names(reformDat) == "value"] <- namesForm
  }else{
    reformDat[ , ((namesForm)) := data.table::tstrsplit(value, ":", fixed = TRUE)][
      , value := NULL]

    if("DP" %in% namesForm){
      reformDat[ , DP := as.numeric(DP)]
      reformDat[is.na(DP), DP := 0]
    }
    if("RO" %in% namesForm){
      reformDat[ , RO := as.numeric(RO)]
      reformDat[is.na(RO), RO := 0]
    }
    if("AO" %in% namesForm){
      reformDat[ , AO := as.numeric(AO)]
      reformDat[is.na(AO), AO := 0]
    }

  }

  # remove not crappy genotype
  # reformDat <- reformDat[!(GT %in% c("./.", ".")), ]
  reformDat <- reformDat[(GT %in% c("./.", ".")), GT := NA]

  # make a  unique marker for each SNP using locus name and SNP position
  # setkey(reformDat, `#CHROM`, POS)
  # implement fread: retained the original column name `#CHROM` vs `#CHROM`
  data.table::setkey(reformDat, `#CHROM`, POS)

  # return the parsed data frame
  return(reformDat)
}
