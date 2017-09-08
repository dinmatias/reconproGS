#' Converts fastsimcoal2 simulated DNA data to a vcf format
#'
#' This parses an arlequin formatted output of fastsimcoal2. Then, it converts
#' the file into a vcf. The fastsimcoal2 simulates alleles, this function
#' assembles those alleles into individual using the given ploidy. Individual
#' genotypes are labeled by individual number and population.
#'
#'
#' @param filePath full path and name of the arp file
#' @param outFile full path and name of the desired outfile file
#' @param ploidy ploidy of the organism. Default is 2 (diploid)
#' @return A data.table with structure similar to a vcf
#' @examples
#' filePath <- "D:/OneDrive/uniFiles/coalescentSimu/snputility/simuRADseq_2Pop/simuRADseq/simuRADseq_1_1.arp"
#' ploidy <- 2
#' x <- arp2vcf(filePath = filePath, ploidy = ploidy)

arp2vcf <- function(filePath = NULL, outFile = NULL, ploidy = 2){

  # from: https://stackoverflow.com/questions/13673894/suppress-nas-in-paste
  paste3 <- function(..., sep =", "){
    L <- list(...)
    L <- lapply(L, function(x) {x[is.na(x)] <- ""; x})
    ret <- gsub(paste0("(^", sep, "|", sep, "$)"), "", gsub(paste0(sep, sep),sep, do.call(paste, c(L,list(sep=sep)))))
    is.na(ret) <- ret==""
    ret
  }

  dataIn <- readLines(filePath)

  stPolyInfo <- grep(x = dataIn, pattern = "#Number")

  numChrom <- as.numeric(tail(unlist(tstrsplit(dataIn[stPolyInfo], " ")),
                              n = 1))

  enPolyInfo <- grep(x = dataIn, pattern = paste("chromosome", numChrom))

  polyInfo <- dataIn[(stPolyInfo + 2):(enPolyInfo + 1)]

  polyInfo <- polyInfo[-grep(x = polyInfo, pattern = "# 0")]

  chromInfo <- polyInfo[seq(1, length(polyInfo), 2)]

  tempInfo <- tstrsplit(chromInfo, " ")


  chromName <- paste(tempInfo[[6]], tempInfo[[7]], sep = "_")

  numSNP <- as.numeric(tempInfo[[2]])


  chromName <- unlist(mapply(FUN = function(z, y) rep(z, y), z = chromName, y = numSNP), use.names = FALSE)

  # rm(tempInfo)
  # rm(numSNP)
  # rm(chromInfo)
  # rm(stPolyInfo)
  # rm(enPolyInfo)
  # rm(numChrom)

  posInfo <-polyInfo[seq(2, length(polyInfo), 2)]

  posInfo <- gsub(x = posInfo, pattern = "#|,", replacement = "")

  posInfo <- unlist(strsplit(posInfo, " "), use.names = FALSE)

  snpData <- data.table::data.table("#CHROM" = chromName, "POS" = posInfo)

  rm(polyInfo)
  rm(posInfo)

  sampInfo <- vector("list", 3)

  sampInfo[[1]] <- grep(x = dataIn, pattern = "SampleData= \\{")

  numPop <- length(sampInfo[[1]])

  sampInfo[[2]] <- sampInfo[[1]] - 1

  sampInfo[[3]] <- as.numeric(gsub(x = dataIn[sampInfo[[2]]],
                                   pattern = "\\D+", replacement = ""))

  tempD <- matrix(nrow = nrow(snpData), ncol = 0)

  for(i in 1:numPop){

    temp <- data.table::tstrsplit(dataIn[(sampInfo[[1]][i] + 1):
                                           (sampInfo[[1]][i] +
                                              sampInfo[[3]][i])],
                                  "\t")[[3]]

    temp <- t(data.table::setDT(data.table::tstrsplit(temp, "")))

    tempD <- cbind(tempD, temp)

  }

  allelesInfo <- apply(tempD, 1, unique)

  maxAlleles <- max(unlist(lapply(allelesInfo, length)))

  allelesInfo <- lapply(allelesInfo, function(x) { tempNA <- rep(NA, maxAlleles); tempNA[1:length(x)] <- x; return(tempNA)})

  allelesInfo <- data.table::data.table(t(setDT(allelesInfo)))

  newGT <- lapply(1:nrow(tempD),
                  FUN = function(x) {
                    match(tempD[x, ], allelesInfo[x, ])
                    })

  tempD <- do.call(rbind, newGT)

  tempD <- tempD - 1

  rm(newGT)

  allelesInfo[ , REF := V1]

  commandMerg <- paste("allelesInfo[ , ALT := paste3(",
                       paste0(names(allelesInfo)[2:(ncol(allelesInfo) - 1)],
                              collapse = ", "),
                       ")]")

  eval(parse(text = commandMerg))

  allelesInfo <- allelesInfo[ , .(REF, ALT)]

  temp <- mapply(FUN = function(x,y) paste(tempD[,x], tempD[,y], sep = "/"),
                 x = seq(1, sum(sampInfo[[3]]), ploidy),
                 y = seq(2, sum(sampInfo[[3]]), ploidy))


  colnames(temp) <- unlist(lapply
                           (1:numPop, function(x){
                             paste("Ind", 1:(sampInfo[[3]][x]/ploidy),
                                   "_Pop", i, sep = "")}))

  snpData[ , ID := "-"]

  snpData <- cbind(snpData, allelesInfo)

  snpData[ , QUAL := "-"][ , FILTER := "PASS"][
    , INFO := "-"][ , FORMAT := "GT"]

  snpData <- cbind(snpData, temp)

  if(!is.null(outFile)){
    data.table::fwrite(testvcf, outFile,
                       quote = FALSE, row.names = FALSE, sep = "\t")

  }

  return(snpData)
}
