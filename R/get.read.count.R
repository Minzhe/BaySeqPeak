###############################################################################
###                          get.read.count.R                               ###
###############################################################################

################  get single bin reads count  ##################
.get.reads.count <- function(i.batch, scan.param, anno, anno.batch.id, sample.bam.info, bam.chrs) {
      
      ### preparation
      anno.batch <- .get.gene.anno(i.batch, anno, anno.batch.id)        # batch annotation
      check.points <- .get.check.points(anno.batch, scan.param$sliding.step)       # location of checkpoint in this batch gene
      num.check.points <- length(check.points)
      num.bams <- length(sample.bam.info$bam)
      
      ### get counts
      read.count <- data.frame(matrix(0, nrow = num.check.points, ncol = num.bams), 
                               batch = rep(i.batch, num.check.points), 
                               check.point = check.points)
      if (sum(bam.chrs == anno.batch$chr) > 0) {
            for (i in 1:num.bams) {
                  # get read count for every bam file
                  read.count[,i] <- .get.check.points.reads.count(sample.bam.info$bam[i], anno.batch, check.points, scan.param)
            }
      } else {
            read.count <- read.count[0,]
      }
      return(read.count)
}


###################  get single batch gene annotation  ########################
.get.gene.anno <- function(i.batch, anno, anno.batch.id) {
      
      ### extract batch annotation
      anno.batch <- anno[anno.batch.id == i.batch, c(1,3:6)]
      anno.unique <- unique(anno.batch)
      
      # extract information
      strand <- as.character(anno.unique[1,4])
      chr <- as.character(anno.unique[1,1])
      left <- min(anno.unique$start)
      right <- max(anno.unique$stop)
      intervals <- anno.unique[,2:3] - left + 1
      gene <- as.character(anno.unique[1,5])
      DNA.length <- right - left + 1
      
      # prepare DNA2RNA
      DNA2RNA <- rep(0, DNA.length)       # initialize with 0 (non-transcribed)
      for (i in 1:length(intervals[,1])) {
            DNA2RNA[intervals[i,1]:intervals[i,2]] <- 1     # mark transcribed region
      }
      exome.length <- sum(DNA2RNA)        # this is actual exome length
      DNA2RNA <- cumsum(DNA2RNA) * DNA2RNA
      
      # summarize result
      batch_anno <- list(gene = gene,
                         chr = chr,
                         strand = strand,
                         left = left,
                         right = right,
                         DNA2RNA = DNA2RNA,
                         DNA.length = DNA.length,
                         exome.length = exome.length)
      
      return(batch_anno)
}



################  get check point for each bin  ##################
.get.check.points <- function(anno.batch, sliding.step){
      
      num.points <- 1 + ceiling(anno.batch$exome.length / sliding.step)       # number of points
      check.points <- as.integer(round(seq(from = 1, to = anno.batch$exome.length, length.out = num.points)))      # generate points
      
      return(check.points)
}



###################  get single bin reads count  #####################
.get.check.points.reads.count<- function(bam, anno.batch, check.points, scan.param){
      
      ### prepare bam parameters
      which <- RangesList(IRanges(anno.batch$left, anno.batch$right))
      names(which) <- anno.batch$chr
      what <- c("strand", "pos", "mapq", "qwidth")
      param <- ScanBamParam(which = which, what = what)
      
      ### read bam file
      report <- scanBam(bam, param = param)
      pos <- report[[1]]$pos - anno.batch$left + 1
      strand <- report[[1]]$strand
      mapq <- report[[1]]$mapq
      read.length <- round(median(report[[1]]$qwidth))
      
      ### fileter mapq and negative pos
      mapq[which(is.na(mapq))] <- 255
      idx <- which(mapq >= scan.param$min.mapq)
      pos <- pos[idx]
      strand <- strand[idx]
      
      ### fileter negative pos
      idx <- which(pos > 0)
      pos <- pos[idx]
      strand <- strand[idx]
      
      ### convert pos into RNA
      RNA.pos <- anno.batch$DNA2RNA[pos]
      idx <- which(((RNA.pos > 0) + !is.na(RNA.pos)) == 2)  # transcribed region
      RNA.pos <- RNA.pos[idx]
      strand <- strand[idx]
      
      ### divide into strand
      pos.pos <- RNA.pos[which(strand == "+")]
      neg.pos <- RNA.pos[which(strand == "-")]
      
      ### shift
      pos.pos <- pos.pos + round(scan.param$fragment.length/2)
      neg.pos <- neg.pos + read.length - round(scan.param$fragment.length/2)
      
      ### merge the two
      pos <- c(pos.pos, neg.pos)
      pos <- pos[pos > 0 & pos < anno.batch$exome.length]
      
      ### get direct count
      check.points.count <- check.points
      pos.table <- table(pos)
      
      ### smooth
      if (scan.param$remove.local.anomalities) pos.table <- .remove.local.anomalities(pos.table)
      
      ### get count
      pos.mapped <- as.numeric(names(pos.table))
      for (i in 1:length(check.points)) { 
            idx <- which(abs(check.points[i] - pos.mapped) * 2 < scan.param$window.width)
            check.points.count[i] <- sum(pos.table[idx])
      }
      
      return(check.points.count)
}


.remove.local.anomalities <- function(pos.table){
      
      # parameters
      max.background.fold.increase <- 4
      background.window <- 200
      
      
      # get position
      pos.mapped <- as.numeric(names(pos.table))
      num.pos.mapped <- length(pos.mapped)
      
      # prepare new table
      new.table <- pos.table
      
      # filter
      for (i in 1:num.pos.mapped) {
            idx <- which(abs(pos.mapped[i] - pos.mapped) < (background.window / 2))
            else.count <- sum(pos.table[idx]) - sum(pos.table[i])
            max.background <- round(else.count * max.background.fold.increase / background.window)
            new.table[i] <- max(1, min(pos.table[i], max.background))
      }
      
      return(new.table)
}
