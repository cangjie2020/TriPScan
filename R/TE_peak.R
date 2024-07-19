#' To calculate the translation efficiency of differentially expressed peaks (diff-peak) genes
#'
#' @param RNA_count Transcriptome result matrix, with the first column as gene ID/transcript ID,
#'    followed by columns for gene counts in different samples.
#' @param CDS_count ribo-seq results matrix, with the first column being gene ID/transcript ID,
#'    followed by several columns representing the RPF of different samples
#' @param gene The output of geneinfo, a Genome annotation matrix
#' @param peak_pic The output of tRNA_diff_peak
#' @param condition Two strings, The names of the experimental group and the control group.
#'    (e.g. condition <- c("Heart","Liver"))
#' @param replicates If the number of biological replicates is the same for both the experimental and control groups,
#'    input only one number. If they are different, input the number of biological replicates for both the
#'    experimental and control groups separately (e.g replicates = c(2,3)).
#' @param gene_id If the first column of RNA_count and CDS_count is gene ID, then gene_id = T.
#'    If it is transcript ID, then gene_id = F.default is gene_id = T
#' @param min_RNA The minimum threshold for the average gene expression level in the transcriptomic matrix.
#' @param min_CDS The minimum threshold for the average RPF level in the ribo-seq results matrix.
#' @param color The color of the bar chart
#'
#'
#' @examples
#'
#' ## RNA<-read.table("../CDS/RNA.txt",sep = "\t",header = T)
#' ## CDS<-read.table("../CDS/CDS.txt",sep = "\t",header = T)
#' ## TE_pic<- TE_peak(RNA_count = RNA,CDS_count = CDS,gene = geneinfo,peak_pic = peak_pic ,condition = condition ,replicates = 2 ,gene_id = F)
#'
#' @import ggplot2
#' @export
#'
TE_peak <- function(RNA_count, CDS_count, gene, peak_pic, condition, replicates,
                    gene_id = T, min_RNA = 10, min_CDS = 10, color = c("#33adff", "#ff4d4d")) {
  output <- list()
  gene <- gene[which(gene$transcript_biotype == "protein_coding"), ]
  gene <- gene[-which(gene$chrom == "MT"), ]

  if (gene_id == F) {
    RNA_count <- dplyr::inner_join(RNA_count, gene[, c(2, 5)], by = "tx_name")
    RNA_count <- RNA_count[rowMeans(RNA_count[, 2:(ncol(RNA_count) - 1)]) >= min_RNA, ]
    RPK <- sweep(RNA_count[, 2:(ncol(RNA_count) - 1)], 1, RNA_count[, ncol(RNA_count)], "/")
    total_RPK <- colSums(RPK)
    TPM <- sweep(RPK, 2, total_RPK, "/") * 1e6
    RNA_count[, 2:(ncol(RNA_count) - 1)] <- TPM

    CDS_count <- dplyr::inner_join(CDS_count, gene[, c(2, 6)], by = "tx_name")
    CDS_count <- CDS_count[rowMeans(CDS_count[, 2:(ncol(CDS_count) - 1)]) >= min_CDS, ]
    RPK <- sweep(CDS_count[, 2:(ncol(CDS_count) - 1)], 1, CDS_count[, ncol(CDS_count)], "/")
    total_RPK <- colSums(RPK)
    TPM <- sweep(RPK, 2, total_RPK, "/") * 1e6
    CDS_count[, 2:(ncol(CDS_count) - 1)] <- TPM
  }

  if (gene_id == T) {
    gene <- gene[order(gene$cds_len, decreasing = T), ]
    gene <- gene[!duplicated(gene$gene_id), ]
    RNA_count <- dplyr::inner_join(RNA_count, gene[, c(1, 2, 5)], by = "gene_id")
    RNA_count <- RNA_count[, -1]
    RNA_count <- RNA_count[, c(ncol(RNA_count) - 1, 1:(ncol(RNA_count) - 2), ncol(RNA_count))]

    RNA_count <- RNA_count[rowMeans(RNA_count[, 2:(ncol(RNA_count) - 1)]) >= min_RNA, ]
    RPK <- sweep(RNA_count[, 2:(ncol(RNA_count) - 1)], 1, RNA_count[, ncol(RNA_count)], "/")
    total_RPK <- colSums(RPK)
    TPM <- sweep(RPK, 2, total_RPK, "/") * 1e6
    RNA_count[, 2:(ncol(RNA_count) - 1)] <- TPM

    CDS_count <- dplyr::inner_join(CDS_count, gene[, c(1, 2, 6)], by = "gene_id")
    CDS_count <- CDS_count[, -1]
    CDS_count <- CDS_count[, c(ncol(CDS_count) - 1, 1:(ncol(CDS_count) - 2), ncol(CDS_count))]

    CDS_count <- CDS_count[rowMeans(CDS_count[, 2:(ncol(CDS_count) - 1)]) >= min_CDS, ]
    RPK <- sweep(CDS_count[, 2:(ncol(CDS_count) - 1)], 1, CDS_count[, ncol(CDS_count)], "/")
    total_RPK <- colSums(RPK)
    TPM <- sweep(RPK, 2, total_RPK, "/") * 1e6
    CDS_count[, 2:(ncol(CDS_count) - 1)] <- TPM
  }
  RNA_CDS <- dplyr::full_join(RNA_count[, 1:(ncol(RNA_count) - 1)], CDS_count[, 1:(ncol(CDS_count) - 1)], by = "tx_name")
  RNA_CDS[is.na(RNA_CDS)] <- 0
  row.names(RNA_CDS) <- RNA_CDS$tx_name
  TE <- RNA_CDS[, ncol(RNA_count):ncol(RNA_CDS)] / RNA_CDS[, 2:(ncol(RNA_count) - 1)]
  TE[] <- sapply(TE, function(x) {
    x[is.nan(x)] <- 0
    x[is.infinite(x)] <- 0
    x
  })
  TE <- TE[rowSums(TE) != 0, ]
  output[["TE_res"]] <- TE

  peak_gene <- peak_pic[[3]]
  TE_pic <- TE[which(row.names(TE) %in% peak_gene$tx_name), ]
  TE_name <- paste("TE", condition[1], "VS", condition[2], sep = "_")
  output[[TE_name]] <- TE_pic
  TE_name1 <- paste(condition[1], "VS", condition[2], sep = "_")
  output[[TE_name1]] <- list()

  for (i in seq_len(nrow(TE_pic))) {
    pic_tmp <- TE_pic[i, ]
    pic_tmp <- as.data.frame(t(pic_tmp))
    if (length(replicates) == 2) {
      mean1 <- mean(pic_tmp[1:replicates[1], ])
      mean2 <- mean(pic_tmp[(replicates[1] + 1):nrow(pic_tmp), ])
      sd1 <- sd(pic_tmp[1:replicates[1], ])
      sd2 <- sd(pic_tmp[(replicates[1] + 1):nrow(pic_tmp), ])
      pic_tmp$Mean <- c(rep(mean1, replicates[1]), rep(mean2, replicates[2]))
      pic_tmp$Sd <- c(rep(sd1, replicates[1]), rep(sd2, replicates[2]))
      pic_tmp$condition <- c(rep(condition[1], replicates[1]), rep(condition[2], replicates[2]))
      t_test_result <- t.test(pic_tmp[1:replicates[1], 1], pic_tmp[(replicates[1] + 1):nrow(pic_tmp), 1], alternative = "two.sided", var.equal = TRUE)
    } else if (length(replicates) == 1) {
      mean1 <- mean(pic_tmp[1:replicates, ])
      mean2 <- mean(pic_tmp[(replicates + 1):nrow(pic_tmp), ])
      sd1 <- sd(pic_tmp[1:replicates, ])
      sd2 <- sd(pic_tmp[(replicates + 1):nrow(pic_tmp), ])
      pic_tmp$Mean <- c(rep(mean1, replicates), rep(mean2, replicates))
      pic_tmp$Sd <- c(rep(sd1, replicates), rep(sd2, replicates))
      pic_tmp$condition <- c(rep(condition[1], replicates), rep(condition[2], replicates))
      t_test_result <- t.test(pic_tmp[1:replicates, 1], pic_tmp[(replicates + 1):nrow(pic_tmp), 1], alternative = "two.sided", var.equal = TRUE)
    }
    text <- paste("P = ", signif(t_test_result$p.value, 3), sep = "")
    colnames(pic_tmp)[1] <- "count"
    P <- ggplot(pic_tmp, aes(x = condition, y = Mean, fill = condition)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "#8c8c8c", size = 0.5) +
      geom_errorbar(aes(ymin = Mean - Sd, ymax = Mean + Sd), position = position_dodge(1), width = .5) +
      geom_point(aes(x = condition, y = count), alpha = 1) +
      ylab("TE") +
      annotate("text", x = 1.35, y = Inf, label = text, vjust = 1.5, hjust = -0.1, size = 7) +
      scale_fill_manual(values = color) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, color = "black"), ) +
      theme(axis.text.y = element_text(size = 20, color = "black"), ) +
      theme(axis.title.x = element_text(size = 0, vjust = 0.5, hjust = 0.5, color = "black")) +
      theme(axis.title.y = element_text(size = 20, vjust = 0.5, hjust = 0.5, color = "black")) +
      theme(panel.grid.major = element_blank()) +
      theme(panel.grid = element_blank()) +
      theme(legend.key = element_blank()) +
      theme(legend.title = element_text(size = 0, face = "bold", color = "black")) +
      theme(legend.text = element_text(size = 20, color = "black"))
    name2 <- gene[which(gene$tx_name == row.names(TE_pic)[i]), 9]
    name2 <- as.character(unique(name2))
    output[[TE_name1]][[name2]] <- P
  }


  peak_gene <- peak_pic[[5]]
  TE_pic <- TE[which(row.names(TE) %in% peak_gene$tx_name), ]
  TE_name <- paste("TE", condition[2], "VS", condition[1], sep = "_")
  output[[TE_name]] <- TE_pic
  TE_name1 <- paste(condition[2], "VS", condition[1], sep = "_")
  output[[TE_name1]] <- list()

  for (i in seq_len(nrow(TE_pic))) {
    pic_tmp <- TE_pic[i, ]
    pic_tmp <- as.data.frame(t(pic_tmp))
    if (length(replicates) == 2) {
      mean1 <- mean(pic_tmp[1:replicates[1], ])
      mean2 <- mean(pic_tmp[(replicates[1] + 1):nrow(pic_tmp), ])
      sd1 <- sd(pic_tmp[1:replicates[1], ])
      sd2 <- sd(pic_tmp[(replicates[1] + 1):nrow(pic_tmp), ])
      pic_tmp$Mean <- c(rep(mean1, replicates[1]), rep(mean2, replicates[2]))
      pic_tmp$Sd <- c(rep(sd1, replicates[1]), rep(sd2, replicates[2]))
      pic_tmp$condition <- c(rep(condition[1], replicates[1]), rep(condition[2], replicates[2]))
      t_test_result <- t.test(pic_tmp[1:replicates[1], 1], pic_tmp[(replicates[1] + 1):nrow(pic_tmp), 1], alternative = "two.sided", var.equal = TRUE)
    } else if (length(replicates) == 1) {
      mean1 <- mean(pic_tmp[1:replicates, ])
      mean2 <- mean(pic_tmp[(replicates + 1):nrow(pic_tmp), ])
      sd1 <- sd(pic_tmp[1:replicates, ])
      sd2 <- sd(pic_tmp[(replicates + 1):nrow(pic_tmp), ])
      pic_tmp$Mean <- c(rep(mean1, replicates), rep(mean2, replicates))
      pic_tmp$Sd <- c(rep(sd1, replicates), rep(sd2, replicates))
      pic_tmp$condition <- c(rep(condition[1], replicates), rep(condition[2], replicates))
      t_test_result <- t.test(pic_tmp[1:replicates, 1], pic_tmp[(replicates + 1):nrow(pic_tmp), 1], alternative = "two.sided", var.equal = TRUE)
    }
    text <- paste("P = ", signif(t_test_result$p.value, 3), sep = "")
    colnames(pic_tmp)[1] <- "count"
    P <- ggplot(pic_tmp, aes(x = condition, y = Mean, fill = condition)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "#8c8c8c", size = 0.5) +
      geom_errorbar(aes(ymin = Mean - Sd, ymax = Mean + Sd), position = position_dodge(1), width = .5) +
      geom_point(aes(x = condition, y = count), alpha = 1) +
      ylab("TE") +
      annotate("text", x = 1.35, y = Inf, label = text, vjust = 1.5, hjust = -0.1, size = 7) +
      scale_fill_manual(values = color) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, color = "black"), ) +
      theme(axis.text.y = element_text(size = 20, color = "black"), ) +
      theme(axis.title.x = element_text(size = 0, vjust = 0.5, hjust = 0.5, color = "black")) +
      theme(axis.title.y = element_text(size = 20, vjust = 0.5, hjust = 0.5, color = "black")) +
      theme(panel.grid.major = element_blank()) +
      theme(panel.grid = element_blank()) +
      theme(legend.key = element_blank()) +
      theme(legend.title = element_text(size = 0, face = "bold", color = "black")) +
      theme(legend.text = element_text(size = 20, color = "black"))
    name2 <- gene[which(gene$tx_name == row.names(TE_pic)[i]), 9]
    name2 <- as.character(unique(name2))
    output[[TE_name1]][[name2]] <- P
  }
  return(output)
}
