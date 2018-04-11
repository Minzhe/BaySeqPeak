######################################################################
###                          parse.bam.R                           ###
######################################################################
library(Rsamtools)

#####  index bam file if.bai not found  #####
.index.bam <- function(ip.bam, input.bam, treated.ip.bam, treated.input.bam) {
      bam.files <- c(ip.bam, input.bam, unlist(treated.ip.bam), unlist(treated.input.bam))
      for (i in 1:length(bam.files)) {
            if (!file.exists(paste(bam.files[i], '.bai', sep=""))) {
                  cat("*** indexing bam file ", bam.files[i], " ...\n", sep = "")
                  indexBam(bam.files[i])
            }
      }
}


#####  get bam file sample id  #####
.get.sample.bam.info <- function(ip.bam, input.bam, treated.ip.bam, treated.input.bam) {
      
      cat("Parsing sample bam file information ... ")
      
      ### number of samples
      num.ip <- length(ip.bam)
      num.input <- length(input.bam)
      num.treated.ip <- length(unlist(treated.ip.bam))
      num.treated.input <- length(unlist(treated.input.bam))
      
      ### generate id: id_untreated_ip, id_untreated_input, id_treated_ip, id_treated_input
      id.untreated.ip <- 1:num.ip
      id.untreated.input <- (num.ip + 1):(num.ip + num.input)
      id.treated.ip <- (num.ip + num.input + 1):(num.ip + num.input + num.treated.ip)
      id.treated.input <- (num.ip + num.input + num.treated.ip + 1):(num.ip + num.input + num.treated.ip + num.treated.input)
      
      ### generate treatment levels
      if (is.list(treated.ip.bam) & is.list(treated.input.bam) & length(treated.ip.bam) > 1 & length(treated.input.bam) > 1) {
            treated.ip.condition.num <- lapply(treated.ip.bam, FUN = length)
            treated.ip.condition <- as.vector(unlist(mapply(FUN = rep, letters[1:length(treated.ip.condition.num)], treated.ip.condition.num)))
            treated.input.condition.num <- lapply(treated.input.bam, FUN = length)
            treated.input.condition <- as.vector(unlist(mapply(FUN = rep, letters[1:length(treated.input.condition.num)], treated.input.condition.num)))
      } else {
            treated.ip.condition <- rep(1, num.treated.ip)
            treated.input.condition <- rep(1, num.treated.input)
      }
      
      ### together id
      id.ip <- c(id.untreated.ip, id.treated.ip)
      id.input <- c(id.untreated.input, id.treated.input)
      
      ### sample names
      sample.tags <- c(rep("untreated.ip", num.ip),
                       rep("untreated.input", num.input),
                       paste("treated.ip.condition", treated.ip.condition, sep = "."),
                       paste("treated.input.condition", treated.input.condition, sep = "."))
      
      ### bam file read length
      bam <- c(ip.bam, input.bam, unlist(treated.ip.bam), unlist(treated.input.bam))
      
      
      ### result
      bam.info <- list(bam = bam,
                       sample.tags = sample.tags,
                       ip = id.ip,
                       input = id.input,
                       untreated.ip = id.untreated.ip,
                       untreated.input = id.untreated.input,
                       treated.ip = id.treated.ip,
                       treated.input = id.treated.input,
                       treated.ip.condition = treated.ip.condition,
                       treated.input.condition = treated.input.condition)
      
      cat("OK\n")
      return(bam.info)
}


#####  get bam file chr info  #####
.get.bam.chrs<- function(file){
      
      param <- ScanBamParam(what = "rname")
      bam <- scanBam(file, param = param)
      chr <- levels(bam[[1]]$rname)       # RNAME in SAM, get reference sequence name
      
      return(chr)
}

