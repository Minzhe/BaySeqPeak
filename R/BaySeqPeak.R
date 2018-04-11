###########################################################################
###                            BaySeqPeak.R                             ###
###########################################################################

BaySeqPeak <- function(ip.bam,
                           input.BAM,
                           treated.ip.bam,
                           treated.input.bam,
                           gtf.path,
                           window.width = 200,
                           sliding.step = 30,
                           fragment.length = 100,
                           min.mapq = 30, 
                           remove.local.anomalities = TRUE) {
      
      ######################   Wrap parameters  #############################
      
      scan.param <- list(window.width = window.width,
                         sliding.step = sliding.step,
                         fragment.length = fragment.length,
                         min.mapq = min.mapq,
                         remove.local.anomalities = TRUE)
      
      
      ########################   prepartion   ###########################
      
      ### read gene annotation
      anno <- .parse.gtf(gtf.path)
      anno.batch.id <- .divide.anno.into.batches(anno)
      
      ### index bam files if .bai file is not found
      .index.bam(ip.bam, input.bam, treated.ip.bam, treated.input.bam)
      
      ### get sample information
      sample.bam.info <- .get.sample.bam.info(ip.bam, input.bam, treated.ip.bam, treated.input.bam)
      bam.chrs <- .get.bam.chrs(ip.bam[1])
      
      cat("Collecting reads count (slow step) ...")
      num.batches <- max(anno.batch.id)
      num.groups <- ceiling(num.batches/170)
      group.bar <- round(seq(from = 1, to = num.batches+1, length.out = num.groups+1))
      
      ### get space
      bin.read.count <- .get.reads.count(1, scan.param, anno, anno.batch.id, sample.bam.info, bam.chrs)
      bin.read.count.format <- bin.read.count <- bin.read.count[0,]
      
      
      ########################   get bin read count   ###########################
      for (i in 1:num.groups){
            cat(as.character(signif(i*100/num.groups, digits = 3)),"%\n", sep = "")
            batches <- group.bar[i]:(group.bar[i+1]-1)
            temp.read.count <- bin.read.count.format
            temp <- mclapply(1:length(batches), FUN = function(x) .get.reads.count(batches[x], scan.param, anno, anno.batch.id, sample.bam.info, bam.chrs))
            temp.read.count <- do.call("rbind", temp)
            bin.read.count <- rbind(bin.read.count, temp.read.count)
      }
      bin.read.count <- data.frame(bin.read.count)
      
      
      ########################   call peaks   ###########################
      
      
      return(list(read.count = bin.read.count, sample.id = sample.bam.info))
}
                    