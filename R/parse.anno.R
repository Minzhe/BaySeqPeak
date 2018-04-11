#######################################################################
###                          parse.anno.R                           ###
#######################################################################
library(GenomicFeatures)

#################  read gtf file  ##################
.parse.gtf <- function(gtf.path){
      
      op.o <- options(); options(warn = -1)
      
      ### make TxDB from GTF file
      cat("Reading GTF annotation file ...\n")
      txdb <- makeTxDbFromGFF(gtf.path, format = "gtf")
      ID <- keys(txdb, "TXID")
      gtf <- select(txdb, ID, c(columns(txdb))[c(7:8,12:14)], "TXID")       # TXID, EXONCHROM, EXONSTRAND, EXONSTART, EXONEND, GENEID
      
      ### make the annotation
      gtf <- cbind(gtf, EXON = NA)
      gtf <- gtf[,c(2,7,4,5,3,6,1)]
      colnames(gtf) <- c("chr", "feature", "start", "stop", "strand", "gene", "transcript")    # EXONCHROM, EXON, EXONSTART, EXONEND, EXONSTRAND, GENEID, TXID
      
      options(op.o)
      
      return(gtf)
}


###################  devide annotation  ########################
.divide.anno.into.batches <- function(anno){
    cat("Divide transcriptome into chr-gene-batch sections ... ")
    
    ### divide into chr-gene
    chr_gene <- paste(anno$chr, anno$gene, anno$strand)
    unique_chr_gene <- unique(chr_gene)
    ID <- match(chr_gene, unique_chr_gene)
    gene_chr_ID <- ID
    
    ### group
    numGroups <- ceiling(max(ID)/100)
    group_bar <- round(seq(from = 1, to = max(ID)+1, length.out = numGroups+1))
    
    for (i in 1:numGroups) {
        # print(i)
        id_selected <- which(((ID >= group_bar[i]) + (ID < group_bar[i+1])) == 2)
        anno_small <- anno[id_selected,]
        ID_small <- ID[id_selected]-group_bar[i]+1
        
        ### return result
        ID_small_new <- .divide.anno.into.batches.small(anno_small, ID_small)
        ID_small_new[ID_small_new > max(ID_small)] <- ID_small_new[ID_small_new > max(ID_small)] - max(ID_small) + max(ID)
        ID_small_new[ID_small_new <= max(ID_small)] <- ID_small_new[ID_small_new <= max(ID_small)] + group_bar[i] - 1
        ID[id_selected] <- ID_small_new
    }
    
    cat("OK\n")
    return(as.integer(ID))
}

.divide.anno.into.batches.small <- function(anno, ID) {
      # divide into chr-gene-batch
      num_batch <- max(ID)
      no_gene_Chr_batch <- num_batch
      for (ibatch in 1:no_gene_Chr_batch) {
            gene_chr_batch <- anno[ID == ibatch, ]
            gene_chr_transcript_unique <- unique(as.character(gene_chr_batch$transcript))
            no_transcript <- length(gene_chr_transcript_unique)
            stop_points <- sort(unique(gene_chr_batch$stop))
            
            if ((no_transcript > 1) & (length(stop_points) > 1)) {
                  # if there are more than 1 transcript annotation
                  # find the breaking points
                  possible_break_points <- stop_points[1:(length(stop_points) - 1)]
                  valid_cut <- possible_break_points * 0
                  
                  # check possible_break_points
                  for (ip in 1:length(possible_break_points)) {
                        left_id <- (gene_chr_batch$stop <= possible_break_points[ip])
                        right_id <- (gene_chr_batch$stop > possible_break_points[ip])
                        batch_transcript <- as.character(gene_chr_batch$transcript)
                        after_break <- length(unique(batch_transcript[left_id])) + length(unique(batch_transcript[right_id]))
                        valid_cut[ip] <- (after_break == no_transcript)
                  }
                  
                  if (sum(valid_cut) > 0) {
                        # cut into pieces
                        valid_cut_points <-
                              possible_break_points[valid_cut > 0]
                        current_ID <- gene_chr_batch$stop * 0 + ibatch
                        for (i in 1:length(current_ID)) {
                              piece_number <- sum(gene_chr_batch$stop[i] > valid_cut_points)
                              if (piece_number > 0) {
                                    current_ID[i] <- piece_number + no_gene_Chr_batch
                              }
                        }
                        
                        # update
                        ID[ID == ibatch] <- current_ID
                        no_gene_Chr_batch <- max(ID)
                  }
            }
      }
      return(ID)
}

