#' Identifying differential tRNA-index
#'
#' @param AC_count The path to the Anticodon_counts_raw.txt file output by mim-tRNA-seq
#'    e.g AC_count <- "Anticodon_counts_raw.txt"
#' @param condition a list,Select the order of the experimental and control groups
#'    in the Anticodon_counts_raw.txt file and list them.
#'    e.g condition <- list(N1 = c(1,2), N2 = c(3,4))
#'    Selecting the experimental group as 1 and 2, and the control group as 3 and 4 from the Anticodon_counts_raw.txt file
#' @param replicates a logical, Does the tRNA-seq data have biological replicates? If yes, set replicates = T,
#'    if no, set replicates = F.
#' @param LFC The log-fold change (LFC) for differential tRNA between two groups, default is 0
#' @param Padj Threshold for significant differential tRNAs, default is 0.05
#' @param color Volcano plot color settings for different tRNA-indexes
#'
#'
#' @examples
#' ## AC_count <- "../Anticodon_counts_raw.txt"
#' ## condition <- list(N1 = c(3,4),N2 = c(5,6))
#' ## replicates = T
#' ## tRNA_res <- diff_AC(AC_count,condition,T)
#'
#' @import ggplot2
#' @export

diff_AC <- function(AC_count, condition, replicates, LFC = 0, Padj = 0.05, color = c("orange", "purple", "grey")) {
  output <- list()

  if (file.exists(AC_count)) {
    AC_count <- read.table(AC_count, sep = "\t", header = T)
    AC_count <- tidyr::separate(data = AC_count, col = Anticodon, into = c("tRNA", "AA", "Anticodon"), sep = "-", remove = F)
    AC_count <- AC_count[!grepl("mito", AC_count$tRNA), ]
    AC_count <- AC_count[, -c(1, ncol(AC_count))]

    AA <- read.table("./inst/extdata/AA.txt", sep = "\t", header = T)
    AC_count <- AC_count[which(AC_count$AA %in% AA$AA), ]
    AC_count <- merge(AA, AC_count[, 2:ncol(AC_count)], by = "Anticodon", all = T)
    AC_count[is.na(AC_count)] <- 0

    ordered_AC <- factor(AC_count$Anticodon, levels = AA$Anticodon, ordered = TRUE)
    AC_count <- AC_count[order(ordered_AC), ]
    row.names(AC_count) <- AC_count$Anticodon
    tRNA_index <- AC_count[, 4:ncol(AC_count)]
    tRNA_index <- tRNA_index[-c(11, 12, 15, 36), ]
    for (n in 4:ncol(AC_count)) {
      s <- NULL
      W <- NULL
      if (is.null(s)) s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)
      p <- 1 - s
      for (i in seq(1, 61, by = 4)) {
        W <- c(
          W,
          p[1] * AC_count[i, n] + p[5] * AC_count[i + 1, n],
          p[2] * AC_count[i + 1, n] + p[6] * AC_count[i, n],
          p[3] * AC_count[i + 2, n] + p[7] * AC_count[i, n],
          p[4] * AC_count[i + 3, n] + p[8] * AC_count[i + 2, n]
        )
      }
      W <- W[-c(11, 12, 15, 36)]
      w <- W / sum(W)
      if (sum(w == 0) > 0) {
        ws <- w[w != 0]
        gm <- exp(sum(log(ws)) / length(ws))
        w[w == 0] <- gm
      }
      tRNA_index[, n - 3] <- w
    }
    tRNA_index$Anticodon <- row.names(tRNA_index)
    tRNA_index <- merge(AA, tRNA_index, by = "Anticodon")

    output[["tRNA_index"]] <- tRNA_index
  } else {
    cat("\n")
    stop(sprintf(
      "incorrect data: \"%s\" not found\n\n",
      paste(data, collapse = ", ")
    ))
  }


  if (replicates == F) {
    if (length(condition$N1) > 1 | length(condition$N2) > 1) {
      cat("\n")
      stop(sprintf("replicates > 1"))
    }

    a <- condition$N1 + 3
    b <- condition$N2 + 3

    A <- tRNA_index[, c(3, a, b)]
    row.names(A) <- A$Codon
    A <- A[, -1]
    A <- A * 1000000

    myfactors <- data.frame(Tissue = c("aexp", "control"), TissueRun = colnames(A))
    myfilt <- NOISeq::filtered.data(A, factor = myfactors$Tissue, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
    mydata <- NOISeq::readData(data = myfilt, factors = myfactors)
    mynoiseq <- NOISeq::noiseq(mydata,
      k = 0.5, norm = "rpkm", factor = "Tissue",
      pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = "no"
    )

    res <- as.data.frame(mynoiseq@results[[1]])
    res$codon <- row.names(res)
    colnames(res)[3] <- "log2FoldChange"
    colnames(res)[5] <- "padj"
    res$padj <- 1 - res$padj

    output[["result"]] <- res
  } else if (replicates == T) {
    if (length(condition$N1) < 2 || length(condition$N2) < 2) {
      cat("\n")
      stop(sprintf("replicates < 2"))
    }

    a <- condition$N1 + 3
    b <- condition$N2 + 3
    A <- tRNA_index[, c(3, b, a)]
    row.names(A) <- A$Codon
    A <- A[, -1]
    A <- A * 1000000
    A[] <- lapply(A, function(x) if (is.numeric(x)) as.integer(round(x)) else x)
    num <- length(a)

    condition <- factor(c(rep("control", num), rep("treat", num)), levels = c("control", "treat"))
    colData <- data.frame(condition, row.names = colnames(A))
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = A,
      colData = colData,
      design = ~condition
    )
    dds$condition <- factor(dds$condition, levels = c("control", "treat"))
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    res <- as.data.frame(res)
    res$codon <- row.names(res)
    output[["result"]] <- res
  }

  res <- res[order(res$padj, decreasing = F), ]
  res$codon <- as.factor(ifelse(res$padj < Padj, res$codon, NA))
  threshold <- as.factor(ifelse(res$padj < Padj, ifelse(res$log2FoldChange > LFC, "Up", "Down"), "Not"))

  P <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
    xlab("LFC") +
    ylab("-log10(FDR)") +
    geom_point(size = 4, alpha = 0.7) +
    geom_vline(aes(xintercept = LFC), lty = 4, col = "black", lwd = 0.8) +
    geom_vline(aes(xintercept = -LFC), lty = 4, col = "black", lwd = 0.8) +
    scale_color_manual(values = color, breaks = c("Up", "Down", "Not")) +
    theme_classic(base_size = 16) +
    ggrepel::geom_text_repel(aes(label = codon), size = 6) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20, color = "black"), ) +
    theme(axis.text.y = element_text(size = 20, color = "black"), ) +
    theme(axis.title.x = element_text(size = 20, vjust = 0.5, hjust = 0.5, color = "black")) +
    theme(axis.title.y = element_text(size = 20, vjust = 0.5, hjust = 0.5, color = "black")) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(panel.border = element_rect(color = "black", linewidth = 1.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(size = 20, color = "black"))
  output[["plot"]] <- P
  return(output)
}
