#' Genome annotation matrix
#'
#' @param gtf.path Corresponding species gtf annotation file
#'
#'
#' @examples
#'
#' ## geneinfo <- extractGeneInfo("../Mus_musculus.GRCm39.110.chr.gtf")
#'
#' @import data.table
#'
#' @export
#'
extractGeneInfo <- function(gtf.path) {
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf.path)
  txlen <- GenomicFeatures::transcriptLengths(txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE)
  x <- rtracklayer::import(gtf.path, "gtf")
  x <- as.data.frame(x)
  scols <- c(
    "gene_id", "gene_name", "gene_biotype", "transcript_id",
    "transcript_biotype", "protein_id", "seqnames", "strand"
  )
  scols <- scols[scols %in% colnames(x)]
  txtype <- x[, scols]
  setDT(txtype)
  txtype <- unique(txtype)
  res <- merge(txlen, txtype, by.x = c("gene_id", "tx_name"), by.y = c("gene_id", "transcript_id"))
  res <- as.data.table(res)
  res <- res[order(res$gene_id, res$tx_name, res$protein_id), ]
  res <- res[!duplicated(res$tx_name), ]
  setnames(res, "seqnames", "chrom")
  return(res)
}
