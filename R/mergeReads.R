#' Reads raw fastq files and concatenates the forward and reverse read
#'
#' Merge the forward and reverse reads of fastq files. This basically concatenates
#' the sequence and the read quality and uses the sequence infromation and other
#' metadata of the forward read. The fastqfile should have the following name format:
#' name.F.fq for the forward and name.R.fq for the reverse.
#'
#' @param filePath full path to the directory containing the forward and reverse reads
#' @param fastqNam name of the fastq file. See Details
#' @param outFile name of the fastq file to be generated
#' @return NULL; saves a fastq file instead.
#' @examples
#' mergeReads(filePath = "rawReads/", fasqNam = "EMSE_5102")

mergeReads <- function(filePath, fastqNam, outFile = NULL, chunkLines = 200000){

  if(is.null(outFile)){
    outFile <- gsub(forward, pattern = ".F.fq", replacement = "con.F.fq")

  }

  numLine <- unname(sapply(c(x,y), FUN = function(d) R.utils::countLines(d)))

  skipList <- seq(0, numLine[1], chunkLines)

  for(i in skipList){
    forReads <- data.table::fread(forward, sep = ".", nrow = chunkLines,
                                  header = FALSE, col.names = "FORWARD",
                                  skip = i)

    revReads <- data.table::fread(reverse, sep = ".", nrow = chunkLines,
                                  header = FALSE, col.names = "REVERSE",
                                  skip = i)

    combiReads <- cbind(forReads, revReads)

    numRows <- nrow(combiReads)

    combiReads[c(seq(2, numRows, 4), seq(4, numRows, 4)) ,
               newCol := paste(FORWARD, REVERSE, sep = "")]

    combiReads[c(seq(1, numRows, 4), seq(3, numRows, 4)) ,
               newCol := FORWARD]

    fwrite(data.table(combiReads[, newCol]), file = outFile,
           col.names = FALSE, append = TRUE)
  }

  print(paste(outFile, "written"))
  return(NULL)


}
