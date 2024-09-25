globalVariables(c("rel_pos", "rel_pos_new"))

#' Create a Whole Metagene Plot
#'
#' This function generates a whole metagene plot based on gene annotation and
#' QC data. It calculates the relative distance of reads from the transcription
#' start site and plots the density of read footprints across different regions
#' of transcripts (5'UTR, CDS, and 3'UTR).
#'
#' @param gene_file A character string specifying the path to the gene annotation file
#' which is from pre_longest_trans_info function. Defaults to `NULL`.
#' @param qc_data A data frame containing QC data with transcript information.
#' which is from load_qc_data function. Defaults to `NULL`.
#' @param scale_region A logical value indicating whether to scale the UTR regions
#' relative to the CDS. Defaults to `FALSE`.
#' @param bins Two length of numerical vector for region scale bins. Default NULL.
#' @param geom_line_list A list of arguments to pass to the `geom_line` function
#' in `ggplot2`. Defaults to an empty list.
#' @param facet_wrap_list A list of arguments to pass to the `facet_wrap` function
#' in `ggplot2`. Defaults to an empty list.
#'
#' @return A ggplot object representing the whole metagene plot.
#' @importFrom stats median
#'
#' @export
whole_metagene_plot <- function(gene_file = NULL,
                                qc_data = NULL,
                                scale_region = FALSE,
                                bins = NULL,
                                geom_line_list = list(),
                                facet_wrap_list = list()){
  # ===================================================================================================
  # load gene annotation
  # ===================================================================================================
  ganao <- read.delim(gene_file,header = F)
  colnames(ganao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                       "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

  ganao <- ganao[,c("trans_id","5UTR_length","CDS_length","3UTR_length")]

  # scale regions
  if(!is.null(bins)){
    utr5_scale <- bins[1]
    utr3_scale <- bins[2]

    scale_region = TRUE
  }else{
    utr5_scale <- stats::median(ganao$`5UTR_length`)/stats::median(ganao$CDS_length)
    utr3_scale <- stats::median(ganao$`3UTR_length`)/stats::median(ganao$CDS_length)
  }

  # ===================================================================================================
  # offset shift
  # ===================================================================================================
  # if(!is.null(offset_shift)){
  #   offset_df <- offset_shift
  #   offset_df <- offset_df %>%
  #     tidyr::separate_longer_delim(cols = c(readLengths,Offsets),delim = ",")
  # }
  # ===================================================================================================
  # calculate relative distance
  # ===================================================================================================
  sp <- unique(qc_data$sample)

  # loop for each sample
  # x = 1
  purrr::map_df(seq_along(sp),function(x){
    tmp <- subset(qc_data,sample == sp[x])

    tmp_anno <- tmp %>%
      dplyr::left_join(y = ganao,by = "trans_id") %>%
      dplyr::mutate(cds_st = `5UTR_length` + 1,
                    cds_sp = `5UTR_length` + CDS_length) %>%
      dplyr::mutate(rel_pos = dplyr::case_when(trans_pos >= 0 & trans_pos < cds_st ~ trans_pos/`5UTR_length`,
                                               trans_pos >= cds_st & trans_pos < cds_sp ~ (trans_pos - cds_st + 1)/CDS_length + 1,
                                               trans_pos >= cds_sp ~ (trans_pos - cds_sp + 1)/`3UTR_length` + 2)) %>%
      dplyr::mutate(rel_pos_new = dplyr::case_when(rel_pos >=0 & rel_pos < 1 ~ scales::rescale(rel_pos,to = c(1 - utr5_scale,1),from = c(0,1)),
                                                   rel_pos >=2 ~ scales::rescale(rel_pos,to = c(2,2 + utr3_scale),from = c(2,3)),
                                                   .default = rel_pos)) %>%
      dplyr::select(sample,group,counts,rel_pos,rel_pos_new)

    return(tmp_anno)
  }) -> relpos_df

  # add conuts to frequency
  df_tmp <- relpos_df[,c("sample","counts","rel_pos","rel_pos_new")]

  sp <- unique(df_tmp$sample)

  # x = 1
  purrr::map_df(seq_along(sp),function(x){
    tmp <- subset(df_tmp,sample == sp[x])
    new_tmp <- data.frame(sample = sp[x],
                          rel_pos = rep(tmp$rel_pos,tmp$counts),
                          rel_pos_new = rep(tmp$rel_pos_new,tmp$counts))
    return(new_tmp)
  }) -> df_plot

  # ===================================================================================================
  # plot
  # ===================================================================================================
  if(scale_region == TRUE){
    breaks <- c((1 - utr5_scale + 1)/2,1.5,(2 + utr3_scale/2))

    line <- do.call(geom_density,modifyList(list(data = df_plot,
                                                 mapping = aes(x = rel_pos_new)),
                                            geom_line_list))
  }else{
    breaks <- c(0.5,1.5,2.5)

    line <- do.call(geom_density,modifyList(list(data = df_plot,
                                                 mapping = aes(x = rel_pos)),
                                            geom_line_list))
  }

  labels <- c("5'UTR","CDS","3'UTR")


  # plot
  ggplot() +
    # geom_density(aes(x = rel_pos_new)) +
    line +
    geom_vline(xintercept = c(1,2),lty = "dashed") +
    theme_bw() +
    scale_x_continuous(breaks = breaks,labels = labels) +
    theme(axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black",face = "bold"),
          strip.text = element_text(color = "black",size = rel(1),face = "bold"),
          axis.text.x = element_text(color = "black",angle = 45,hjust = 1),
          panel.grid = element_blank()) +
    do.call(facet_wrap,
            modifyList(list(facets = vars(sample),scales = 'fixed'),
                       facet_wrap_list)) +
    ylab("Relative footprint density (AU)") +
    xlab("Normalized transcript length")
}


