get.table.peak.result <- function(PEAK, ANNOTATION, READS_COUNT, ANNOTATION_BATCH_ID, LOCI2PEAK, PPI, FC) {
      
      # initializa the peak reporting
      no_peak <- length(LOCI2PEAK[,1])
      peak_report <- data.frame()
      if (no_peak == 0) {
            return(peak_report)
      } else {
            # get peak
            for (i in 1:no_peak) {
                  if (i %% 100 == 1) print(i)
                  peak_row_id <- LOCI2PEAK[i,]
                  rna_peak <- READS_COUNT$check_points[peak_row_id]
                  
                  # batch id
                  batch_id <- unique(READS_COUNT$batch_id[peak_row_id])
                  ppi <- max(PPI[peak_row_id[1]:peak_row_id[2]])
                  count.table <- READS_COUNT[peak_row_id[1]:peak_row_id[2],1:6]
                  fold_change_bin <- max(apply(count.table, 1, FUN = function(x) sum(x[1:3]+1)/sum(x[4:6]+1)))
                  fold_change <- max(FC[peak_row_id[1]:peak_row_id[2]])
                  # get sig digits
                  ppi <- signif(ppi, digits = 3)
                  # get annotation
                  anno <- .get.gene.anno(batch_id,ANNOTATION,ANNOTATION_BATCH_ID)
                  # get blocks
                  block <- .get.block.from.rna.peak(rna_peak,anno)
                  # save result
                  xls <- data.frame(chr=anno$chr,
                                    chromStart=block$chromStart,
                                    chromEnd=block$chromEnd,
                                    name=anno$gene,
                                    score=signif(1-ppi,digits=2),
                                    strand=anno$strand,
                                    thickStart=block$thickStart,
                                    thickEnd=block$thickEnd,
                                    itemRgb=0,
                                    blockCount=block$blockCount,
                                    blockSizes=block$blockSizes,
                                    blockStarts=block$blockStarts,
                                    ppi=ppi,
                                    fold_change_bin=fold_change_bin,
                                    fold_change=fold_change)
                  # append result
                  peak_report=rbind(peak_report,xls)
            }
            return(peak_report)}
}

get.peak.from.loci <- function(READS_COUNT, ID, MINIMAL_PEAK_LENGTH = 50){
      # get no_line
      no_line <- length(ID)
      
      # start id
      start_id <- which((ID[2:no_line]-ID[1:no_line-1]==1) | ((READS_COUNT$batch_id[2:no_line]!=READS_COUNT$batch_id[1:no_line-1]) & (ID[2:no_line] == TRUE)))
      start_id=start_id + 1
      if (ID[1]==TRUE) { start_id <- c(1,start_id) }
      
      # end id 
      end_id = which((ID[1:no_line-1]-ID[2:no_line] == 1) | ((READS_COUNT$batch_id[1:no_line-1]!=READS_COUNT$batch_id[2:no_line]) & (ID[1:no_line-1] == TRUE)) )
      if (ID[no_line] == TRUE) {end_id <- c(end_id,no_line)}  
      
      # label peaks
      PEAK_LENGTH <- READS_COUNT$check_points[end_id]-READS_COUNT$check_points[start_id]
      good_peak_id <- which(PEAK_LENGTH > MINIMAL_PEAK_LENGTH)
      start_id <- start_id[good_peak_id]
      end_id <- end_id[good_peak_id]
      
      # result
      peak <- cbind(start_id,end_id)
      return(peak)
}

get.peak <- function(READS_COUNT, PPI, CUTOFF = 0.5, MIN_COUNT = 10, MINIMAL_PEAK_LENGTH = 50) {
      PEAK <- list()
      PEAK$ppi <- PPI
      PEAK$Merged <- PPI >= CUTOFF & rowSums(READS_COUNT[,1:6]) > MIN_COUNT
      PEAK$loci2peak_merged <- get.peak.from.loci(READS_COUNT, ID = PEAK$Merged, MINIMAL_PEAK_LENGTH)
      return(PEAK)
}

write.bed <- function(TOTAL_PEAK_RESULT, path) {
      temp <- TOTAL_PEAK_RESULT[,1:12]
      # temp <- apply(temp, 2, as.character)
      names(temp)[1]=paste("#",names(temp)[1])
      write.table(temp, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}