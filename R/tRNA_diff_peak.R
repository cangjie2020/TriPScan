#' Combine the results of EPA_res and diff_AC to identify the differential peaks caused by tRNA.
#'
#' @param EPA_path a list, List the results of EPA_res in the order of the experimental group and the control group.
#'    e.g if replicates = F  EPA_path <- list(S1 = list(EPA_result[[1]]), S2 = list(EPA_result[[2]]))
#'    e.g if replicates = T  EPA_path <- list(S1 = list(EPA_result[[1]],EPA_result[[2]]),S2 = list(EPA_result[[3]],EPA_result[[4]]))
#' @param condition Two strings, The names of the experimental group and the control group.
#'    (e.g. condition <- c("Heart","Liver"))
#' @param replicates a logical, Does the Ribo-seq data have biological replicates? If yes, set replicates = T,
#'    if no, set replicates = F.
#' @param tRNA  The output of diff_AC.
#' @param geneinfo The output of geneinfo, a Genome annotation matrix
#' @param min_diff Set the minimum peak ratio difference for the same gene at the same locus between two group.default is 3
#' @param diff_tRNA_num Set the significant difference in the number of tRNAs. If diff_tRNA_num = "all",
#'    then input all significant different tRNAs. If the value of diff_tRNA_num is greater than
#'    the number of significant different tRNAs, input all significant different tRNAs.
#'    If the value of diff_tRNA_num is less than the number of significant different tRNAs,
#'    input the diff_tRNA_num most significant tRNAs. default is 8
#' @param padj Threshold for significant differential tRNAs, default is 0.05
#' @param LFC The log-fold change (LFC) for differential tRNA between two groups, default is 0
#' @param min_count The minimum number of ribosome-protected fragments (RPFs) per gene, default is 10
#' @param Non_translation_Elongation Exclude peak judgments associated with translation initiation and termination sites.
#'    Exclude peak results from regions downstream of start codons and upstream of stop codons by n codons. default is 5
#' @param peak_method Calculate peak values using two methods: mean, which computes the ratio of reads at each site to
#'    the average reads per site across the gene; and median, which computes the ratio of reads at each site to the
#'    median reads per site across the gene. default is mean
#' @param only_both Whether differential peak sites have reads in both groups, default is F
#' @param bar_color The colors of the bar chart, default is '#1F78B4'
#' @param line_color The colors of differential peak sites, default is '#FE4A01'
#'
#'
#' @examples
#' ## peak_pic <- tRNA_diff_peak(EPA_path = EPA_path,tRNA = tRNA_res,condition = condition,replicates = T,
#' ##            geneinfo = geneinfo,min_diff = 8,diff_tRNA_num = 5)
#'
#' @import data.table
#' @import ggplot2
#'
#' @export
#'
tRNA_diff_peak <- function(EPA_path, condition, replicates, tRNA, geneinfo,
                           min_diff = 3, diff_tRNA_num = 8, padj = 0.05, LFC = 0, min_count = 10,
                           Non_translation_Elongation = 5, peak_method = "mean", only_both = F,
                           bar_color = "#1F78B4", line_color = "#FE4A01") {
  tRNA <- tRNA$result
  gene <- geneinfo
  gene <- gene[which(gene$transcript_biotype == "protein_coding"), ]
  gene <- gene[which(gene$cds_len %% 3 == 0), ]

  output <- list()
  if (replicates == T) {
    if (length(EPA_path$S1) < 2 || length(EPA_path$S2) < 2) {
      cat("\n")
      stop(sprintf("replicates < 2 \n\n"))
    }
    for (n in 1:2) {
      name_list <- c()
      for (i in seq_len(length(EPA_path[[n]]))) {
        A <- EPA_path[[n]][[i]]
        A <- dplyr::inner_join(A, gene[, c(2, 6, 7)], by = "tx_name")
        A$dis <- as.numeric(A$trans.site) - A$utr5_len - 1
        A$dis <- A$dis / 3
        A$dis2 <- A$utr5_len + A$cds_len - 2 - as.numeric(A$trans.site)
        A$dis2 <- A$dis2 / 3
        A <- A[which(A$dis > Non_translation_Elongation & A$dis2 > Non_translation_Elongation), ]

        A <- A[, c("tx_name", "dis", "A_Codon")]
        A1 <- as.data.frame(table(A[, c("tx_name", "dis")]))
        A1 <- A1[-which(A1$Freq == 0), ]
        A1$name <- paste(A1[, 1], A1[, 2], sep = "_")
        A <- unique(as.data.table(A))
        A$dis <- as.numeric(A$dis)
        A1$dis <- as.numeric(as.character(A1$dis))
        A1 <- dplyr::inner_join(A1, A, by = c("tx_name", "dis"))
        assign(paste(condition[n], i, "peak", sep = "_"), A1)
        name_list <- c(name_list, paste(condition[n], i, "peak", sep = "_"))
      }

      peak_codon <- rbind(get(name_list[1]), get(name_list[2]))
      if (length(name_list) >= 3) {
        for (n in 3:length(name_list)) {
          peak_codon <- rbind(peak_codon, get(name_list[n]))
        }
      }
      peak_codon <- peak_codon[, c(1, 2, 4, 5)]
      peak_codon <- unique(as.data.table(peak_codon))

      A1 <- dplyr::full_join(get(name_list[1])[, c(4, 3)], get(name_list[2])[, c(4, 3)], by = "name")
      if (length(name_list) >= 3) {
        for (n in 3:length(name_list)) {
          A1 <- dplyr::full_join(A1, get(name_list[n])[, c(4, 3)], by = "name")
        }
      }
      A1[is.na(A1)] <- 0
      A1$count <- apply(A1[, 2:ncol(A1)], 1, mean)
      A1 <- dplyr::inner_join(A1, peak_codon, by = "name")

      if (peak_method == "mean") {
        c1 <- aggregate(A1$count, by = list(type = A1$tx_name), mean)
      } else if (peak_method == "median") {
        c1 <- aggregate(A1$count, by = list(type = A1$tx_name), median)
      }

      colnames(c1)[1] <- "tx_name"
      A1 <- dplyr::inner_join(A1, c1, by = "tx_name")
      A1$ratio <- A1$count / A1$x
      A1 <- A1[, c("name", "count", "ratio", "A_Codon")]
      A1 <- unique(as.data.table(A1))
      if (n == 1) {
        S1 <- A1
      } else if (n == 2) {
        S2 <- A1
      }
    }
  } else if (replicates == F) {
    if (length(EPA_path$S1) >= 2 || length(EPA_path$S2) >= 2) {
      cat("\n")
      stop(sprintf("replicates > 1 \n\n"))
    }
    for (i in 1:2) {
      A <- EPA_path[[i]][[1]]
      A <- dplyr::inner_join(A, gene[, c(2, 6, 7)], by = "tx_name")
      A$dis <- as.numeric(A$trans.site) - A$utr5_len - 1
      A$dis <- A$dis / 3
      A$dis2 <- A$utr5_len + A$cds_len - 2 - as.numeric(A$trans.site)
      A$dis2 <- A$dis2 / 3
      A <- A[which(A$dis > Non_translation_Elongation & A$dis2 > Non_translation_Elongation), ]

      A <- A[, c("tx_name", "dis", "A_Codon")]
      A1 <- as.data.frame(table(A[, c("tx_name", "dis")]))
      A1 <- A1[-which(A1$Freq == 0), ]
      A1$name <- paste(A1[, 1], A1[, 2], sep = "_")
      A <- unique(as.data.table(A))
      A$dis <- as.numeric(A$dis)
      A1$dis <- as.numeric(as.character(A1$dis))
      A1 <- dplyr::inner_join(A1, A, by = c("tx_name", "dis"))
      if (peak_method == "mean") {
        c1 <- aggregate(A1$Freq, by = list(type = A1$tx_name), mean)
      } else if (peak_method == "median") {
        c1 <- aggregate(A1$Freq, by = list(type = A1$tx_name), median)
      }

      colnames(c1)[1] <- "tx_name"
      A1 <- dplyr::inner_join(A1, c1, by = "tx_name")
      A1$ratio <- A1$Freq / A1$x
      colnames(A1)[3] <- "count"
      A1 <- A1[, c("name", "count", "ratio", "A_Codon")]
      A1 <- unique(as.data.table(A1))
      if (i == 1) {
        S1 <- A1
      } else if (i == 2) {
        S2 <- A1
      }
    }
  }

  count_tmp <- S1
  count_tmp <- as.data.table(count_tmp)[, c("tx_name", "site") := tstrsplit(name, "_")]
  c1 <- aggregate(count_tmp$count, by = list(type = count_tmp$tx_name), sum)
  c1 <- c1[which(c1$x > min_count), ]

  count_tmp <- S2
  count_tmp <- as.data.table(count_tmp)[, c("tx_name", "site") := tstrsplit(name, "_")]
  c2 <- aggregate(count_tmp$count, by = list(type = count_tmp$tx_name), sum)
  c2 <- c2[which(c2$x > min_count), ]
  C <- dplyr::inner_join(c1, c2, by = "type")

  A <- rbind(S1[, c(1, 4)], S2[, c(1, 4)])
  A <- unique(as.data.table(A))
  A1 <- dplyr::full_join(S1[, c(1, 3)], S2[, c(1, 3)], by = "name")
  A1[is.na(A1)] <- 0
  A1$diff <- A1[, 2] - A1[, 3]
  A1 <- as.data.table(A1)[, c("tx_name", "site") := tstrsplit(name, "_")]
  A1 <- A1[which(A1$tx_name %in% C$type), ]

  A1$site <- as.numeric(A1$site)
  A <- dplyr::inner_join(A1, A, by = "name")
  colnames(A)[2] <- paste("ratio", condition[1], sep = "_")
  colnames(A)[3] <- paste("ratio", condition[2], sep = "_")
  output[["diff_peak"]] <- A

  if (only_both == T) {
    A <- A[-which(A$ratio_Heart == 0 | A$ratio_Liver == 0), ]
  } else if (only_both == F) {
    A <- A
  }

  S1$type <- condition[1]
  S2$type <- condition[2]
  pic <- rbind(S1, S2)
  pic <- as.data.table(pic)[, c("tx_name", "site") := tstrsplit(name, "_")]
  pic$site <- as.numeric(pic$site)

  high_tmp <- A[which(A$diff > min_diff), ]
  low_tRNA <- tRNA[which(tRNA$padj < padj & tRNA$log2FoldChange < -LFC), ]
  low_tRNA <- low_tRNA[order(low_tRNA$log2FoldChange), ]

  if (diff_tRNA_num >= nrow(low_tRNA)) {
    high_tmp <- high_tmp[which(high_tmp$A_Codon %in% low_tRNA$codon), ]
  } else if (diff_tRNA_num < nrow(low_tRNA)) {
    high_tmp <- high_tmp[which(high_tmp$A_Codon %in% low_tRNA$codon[1:diff_tRNA_num]), ]
  } else if (diff_tRNA_num == "all") {
    high_tmp <- high_tmp[which(high_tmp$A_Codon %in% low_tRNA$codon), ]
  }

  high_tmp <- high_tmp[order(high_tmp$diff, decreasing = T), ]

  pic1 <- paste(condition[1], "VS", condition[2], sep = "_")
  output[[pic1]] <- list()
  high_tmp <- as.data.frame(high_tmp)
  pic11 <- paste("res", condition[1], "VS", condition[2], sep = "_")
  output[[pic11]] <- high_tmp
  for (i in high_tmp$tx_name) {
    pic_tmp <- pic[which(pic$tx_name == i), ]
    site1 <- as.numeric(high_tmp[which(high_tmp$tx_name == i), 6])
    pic_tmp$A_Codon <- as.factor(ifelse(pic_tmp$site %in% site1, pic_tmp$A_Codon, NA))
    codon1 <- pic_tmp$A_Codon
    codon1 <- na.omit(codon1)
    codon1 <- as.character(codon1)
    codon1 <- unique(codon1)
    P <- ggplot(pic_tmp, aes(x = site, y = count)) +
      geom_vline(xintercept = site1, linetype = "dashed", col = line_color) +
      geom_bar(stat = "identity", position = "dodge", color = bar_color, fill = bar_color) +
      ggrepel::geom_text_repel(aes(label = A_Codon, y = count), size = 6) +
      xlab("Distance from the start codon") +
      ylab("RPF") +
      geom_hline(yintercept = -0.1, col = "black", lwd = 0.2, show.legend = FALSE) +
      theme_bw() +
      theme_classic() +
      theme(axis.title.x=element_text(vjust=2, size=20))+
      theme(axis.title.y=element_text(vjust=2, size=20))+
      theme(axis.text.x = element_text(size = 20, colour = "black"), )  +
      theme(axis.text.y = element_text(size = 20, colour = "black"), ) +
      facet_grid(type ~ ., scales = "free_y") +
      theme(strip.text.y = element_text(size = 20, colour = "black"))
    name2 <- gene[which(gene$tx_name == i), 9]
    name2 <- unique(name2)
    plot_name <- paste(name2, codon1, sep = "_")
    if (length(plot_name) > 1) {
      for (plots in plot_name) {
        output[[pic1]][[plots]] <- P
      }
    } else {
      output[[pic1]][[plot_name]] <- P
    }
  }

  low_tmp <- A[which(A$diff < -min_diff), ]
  high_tRNA <- tRNA[which(tRNA$padj < padj & tRNA$log2FoldChange > LFC), ]
  high_tRNA <- high_tRNA[order(high_tRNA$log2FoldChange), ]

  if (diff_tRNA_num >= nrow(high_tRNA)) {
    low_tmp <- low_tmp[which(low_tmp$A_Codon %in% high_tRNA$codon), ]
  } else if (diff_tRNA_num < nrow(high_tRNA)) {
    low_tmp <- low_tmp[which(low_tmp$A_Codon %in% high_tRNA$codon[1:diff_tRNA_num]), ]
  } else if (diff_tRNA_num == "all") {
    low_tmp <- low_tmp[which(low_tmp$A_Codon %in% high_tRNA$codon), ]
  }

  low_tmp <- low_tmp[order(low_tmp$diff, decreasing = F), ]
  pic2 <- paste(condition[2], "VS", condition[1], sep = "_")
  output[[pic2]] <- list()
  low_tmp <- as.data.frame(low_tmp)
  pic22 <- paste("res", condition[2], "VS", condition[1], sep = "_")
  output[[pic22]] <- low_tmp

  for (i in low_tmp$tx_name) {
    pic_tmp <- pic[which(pic$tx_name == i), ]

    site1 <- as.numeric(low_tmp[which(low_tmp$tx_name == i), 6])
    pic_tmp$A_Codon <- as.factor(ifelse(pic_tmp$site %in% site1, pic_tmp$A_Codon, NA))
    codon1 <- pic_tmp$A_Codon
    codon1 <- na.omit(codon1)
    codon1 <- as.character(codon1)
    codon1 <- unique(codon1)
    P <- ggplot(pic_tmp, aes(x = site, y = count)) +
      geom_vline(xintercept = site1, linetype = "dashed", col = line_color) +
      geom_bar(stat = "identity", position = "dodge", color = bar_color, fill = bar_color) +
      ggrepel::geom_text_repel(aes(label = A_Codon, y = count), size = 6) +
      xlab("Distance from the start codon") +
      ylab("RPF") +
      geom_hline(yintercept = -0.1, col = "black", lwd = 0.2, show.legend = FALSE) +
      theme_bw() +
      theme_classic() +
      theme(axis.title.x=element_text(vjust=2, size=20))+
      theme(axis.title.y=element_text(vjust=2, size=20))+
      theme(axis.text.x = element_text(size = 20, colour = "black"), ) +
      theme(axis.text.y = element_text(size = 20, colour = "black"), ) +
      facet_grid(type ~ ., scales = "free_y") +
      theme(strip.text.y = element_text(size = 20, colour = "black"))
    name2 <- gene[which(gene$tx_name == i), 9]
    name2 <- unique(name2)
    plot_name <- paste(name2, codon1, sep = "_")
    if (length(plot_name) > 1) {
      for (plots in plot_name) {
        output[[pic2]][[plots]] <- P
      }
    } else {
      output[[pic2]][[plot_name]] <- P
    }
  }
  return(output)
}
