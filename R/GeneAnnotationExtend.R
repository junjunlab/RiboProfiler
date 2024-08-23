#' gene_anno_extend
#'
#' @param longest_trans_file A character string specifying the path to the file
#' containing the longest transcript information.
#' @param upstream_scale An integer specifying the scale for extending the 5' UTR
#' region upstream. If not provided, the 5' UTR length is used.
#' @param upstream_extend An integer specifying the number of nucleotides to extend
#' the 5' UTR region upstream. If not provided, the 5' UTR length is used.
#' @param downstream_scale An integer specifying the scale for extending the 3' UTR
#' region downstream. If not provided, the 3' UTR length is used.
#' @param downstream_extend An integer specifying the number of nucleotides to extend
#' the 3' UTR region downstream. If not provided, the 3' UTR length is used.
#' @param output_file A character string specifying the path to the output file.
#' If not provided, the function returns the extended data frame.
#'
#' @return A data frame with extended gene annotation information.
#'
#' @export
gene_anno_extend <- function(longest_trans_file = NULL,
                             upstream_scale = NULL,
                             upstream_extend = NULL,
                             downstream_scale = NULL,
                             downstream_extend = NULL,
                             output_file = NULL){
  # extend region 50nt upstream and downstream
  df <- read.delim(longest_trans_file,header = F)

  # add colnames
  colnames(df) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                    "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

  # add gene start and stop position
  split_rg <- strsplit(df$exon_rg,split = ":")

  gstart <- purrr::map(split_rg,.f = function(x){x[[1]]}) %>% unlist() %>% as.integer()
  gstop <- purrr::map(split_rg,.f = function(x){x[[length(x)]]}) %>% unlist() %>% as.integer()


  # loop to extend exon region
  purrr::map_df(seq_along(split_rg),function(x){
    tmp_anno <- df[x,]
    utr5 <- tmp_anno$`5UTR_length`
    utr3 <- tmp_anno$`3UTR_length`

    # ==========================================================================================
    # deal with upstream
    # ==========================================================================================
    if(!is.null(upstream_scale) & is.null(upstream_extend)){
      tmp_anno$`5UTR_length` <- upstream_scale
      exd <- upstream_scale - utr5
      if(exd < 0){
        message(paste(tmp_anno$gene_name,"UTR5 length biger than",upstream_scale,"!",sep = " "))
      }
      gstart_new <- gstart[x] - exd
    }else if(is.null(upstream_scale) & !is.null(upstream_extend)){
      tmp_anno$`5UTR_length` <- utr5 + upstream_extend
      gstart_new <- gstart[x] - upstream_extend
    }else{
      gstart_new <- gstart[x]
    }

    # ==========================================================================================
    # deal with downstream
    # ==========================================================================================
    if(!is.null(downstream_scale) & is.null(downstream_extend)){
      tmp_anno$`3UTR_length` <- downstream_scale
      exd <- downstream_scale - utr3
      if(exd < 0){
        message(paste(tmp_anno$gene_name,"UTR3 length biger than",downstream_scale,"!",sep = " "))
      }
      gstop_new <- gstop[x] + exd
    }else if(is.null(downstream_scale) & !is.null(downstream_extend)){
      tmp_anno$`3UTR_length` <- utr3 + downstream_extend
      gstop_new <- gstop[x] + downstream_extend
    }else{
      gstop_new <- gstop[x]
    }

    # ==========================================================================================
    # connect exon regions
    # ==========================================================================================
    tmp <- split_rg[[x]]

    if(length(tmp) == 2){
      exon_rg_new <- paste(as.integer(gstart_new),as.integer(gstop_new),sep = ":")
    }else{
      exon_rg_new <- base::paste(as.integer(gstart_new),
                                 base::paste0(tmp[2:(length(tmp) - 1)],collapse = ":"),
                                 as.integer(gstop_new),sep = ":")
    }

    # output
    tmp_anno$exon_rg <- exon_rg_new

    return(tmp_anno)
  }) -> extended_df


  # output
  if(!is.null(output_file)){
    write.table(extended_df,file = output_file,
                col.names = F,quote = F,sep = "\t",row.names = F)
  }

  return(extended_df)
}



globalVariables(c("3UTR_length", "5UTR_length", "CDS_length"))
#' gene_anno2Ribominer
#'
#' @param longest_trans_file A character string specifying the path to the file
#' containing the longest transcript information.
#' @param output_file A character string specifying the path to the output file.
#' If not provided, the function returns the processed data frame.
#'
#' @return A data frame with processed gene annotation information in the format
#' required by RiboMiner.
#'
#' @export
gene_anno2Ribominer <- function(longest_trans_file = NULL,
                                output_file = NULL){
  gene_anno <- read.delim(longest_trans_file,header = F)

  # add colnames
  colnames(gene_anno) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                           "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

  gene_anno <- gene_anno %>%
    dplyr::mutate(CDS_start = as.integer(`5UTR_length` + 1),
                  CDS_stop = as.integer(`5UTR_length` + CDS_length),
                  transcript_length = `5UTR_length` + CDS_length + `3UTR_length`)

  # add gene start and stop position
  split_rg <- strsplit(gene_anno$exon_rg,split = ":")

  gstart <- purrr::map(split_rg,.f = function(x){x[[1]]}) %>% unlist() %>% as.integer()
  gstop <- purrr::map(split_rg,.f = function(x){x[[length(x)]]}) %>% unlist() %>% as.integer()

  gene_anno$gene_start <- gstart
  gene_anno$gene_stop <- gstop

  # re-arrange
  gene_anno_new <- gene_anno[,c("chrom","trans_id",	"strand",	"gene_id",	"gene_name",
                                "gene_start",	"gene_stop",
                                "CDS_start",	"CDS_stop",
                                "CDS_length",	"5UTR_length",	"3UTR_length",
                                "transcript_length")] |>
    dplyr::mutate(transcript_biotype = "protein_coding",.before = "gene_start")

  # output
  write.table(gene_anno_new,file = output_file,
              col.names = T,quote = F,sep = "\t",row.names = F)
}
