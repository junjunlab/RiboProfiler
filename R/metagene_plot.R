globalVariables(c("avergae_exp", "group_name", "mean_exp","framesp","framest"))

#' Metagene Plot Generation
#'
#' This function generates a metagene plot based on normalized expression data across transcripts.
#' It allows for flexible filtering, normalization, and plotting of gene expression data relative
#' to translation start (st) or stop (sp) sites. Data can be plotted at nucleotide or codon resolution.
#'
#' @param longest_trans_file A file path to the longest transcript data, containing gene annotations.
#'        Default is NULL.
#' @param normed_file A file path to the normalized gene expression data. Default is NULL.
#' @param min_counts The minimum count threshold to include a transcript in the analysis.
#'        Default is 64.
#' @param min_cds_length The minimum CDS length required to include a transcript. Default is 600.
#' @param relative_distance A numeric vector specifying the relative distance range from the start
#'        or stop codon to plot. Default is c(-25, 500).
#' @param type Character vector, either "st" for start codon or "sp" for stop codon, to specify the
#'        reference point for relative distances. Default is c("st", "sp") and will be matched in order.
#' @param mode Character vector, either "nt" for nucleotide or "codon" for codon resolution in the plot.
#'        Default is c("nt", "codon") and will be matched in order.
#' @param collapse Logical, whether to collapse data across samples or not. Default is FALSE.
#' @param frame Whether add in-frame information for metagene plot. Default is FALSE.
#' @param geom_line_params List of additional parameters for geom_line, used to customize the line elements
#'        of the plot. Default is an empty list().
#' @param facet_wrap_params List of additional parameters for facet_wrap, used to customize the faceting
#'        of the plot. Default is an empty list().
#'
#' @return A ggplot object representing the metagene analysis plot.
#' @import ggplot2
#' @import dplyr
#' @import purrr
#' @importFrom rlang sym syms
#'
#' @examples
#' \dontrun{# Assuming `longest_trans_file` and `normed_file` are paths to your data files:
#' metagene_plot(
#'   longest_trans_file = "path/to/longest_transcripts.tsv",
#'   normed_file = "path/to/normed_expression.tsv",
#'   min_counts = 64,
#'   min_cds_length = 600,
#'   relative_distance = c(-25, 500),
#'   type = "st",
#'   mode = "nt"
#' )}
#'
#' @export
metagene_plot <- function(longest_trans_file = NULL,
                          normed_file = NULL,
                          min_counts = 64,
                          min_cds_length = 600,
                          relative_distance = c(-25,500),
                          type = c("st","sp"),
                          mode = c("nt","codon"),
                          collapse = FALSE,
                          frame = FALSE,
                          geom_line_params = list(),
                          facet_wrap_params = list()){
  # check args
  mode <- match.arg(mode,c("nt","codon"))
  type <- match.arg(type,c("st","sp"))

  # check whether add frame info
  if(type == "st"){
    if(frame == TRUE){
      vars_f <- rlang::syms(c("relst","framest"))
    }else{
      vars_f <- rlang::sym("relst")
    }

  }else{
    if(frame == TRUE){
      vars_f <- rlang::syms(c("relsp","framesp"))
    }else{
      vars_f <- rlang::sym("relsp")
    }
  }
  # ====================================================================================
  # gene annotation
  # ====================================================================================
  gene_anao <- read.delim(longest_trans_file,header = F)
  colnames(gene_anao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                           "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

  # filter cds length
  gene_lenth <- gene_anao %>%
    dplyr::select(trans_id,`5UTR_length`,CDS_length,`3UTR_length`) %>%
    dplyr::filter(CDS_length >= min_cds_length)

  # filter counts in CDS region
  total_exp <- normed_file %>%
    dplyr::filter(trans_id %in% gene_lenth$trans_id) %>%
    dplyr::left_join(y = gene_lenth,by = "trans_id") %>%
    dplyr::filter(trans_pos >= `5UTR_length` & trans_pos <= `5UTR_length` + CDS_length) %>%
    dplyr::group_by(trans_id) %>%
    dplyr::summarise(total_counts = sum(counts)) %>%
    dplyr::filter(total_counts >= min_counts)

  # ====================================================================================
  # loop deal with for each samples
  # ====================================================================================
  sp <- unique(normed_file$sample)

  # x = 1
  purrr::map_df(seq_along(sp),function(x){
    tmp <- subset(normed_file, sample == sp[x]) %>%
      dplyr::filter(trans_id %in% total_exp$trans_id)

    trans_numbers <- length(unique(tmp$trans_id))

    if(frame == TRUE){
      tmp_1 <- tmp %>%
        dplyr::group_by(sample,group,!!vars_f[[1]],!!vars_f[[2]]) %>%
        dplyr::summarise(avergae_exp = sum(norm_exp)/trans_numbers) %>%
        dplyr::filter(!!vars_f[[1]] >= relative_distance[1] & !!vars_f[[1]] <= relative_distance[2])

      # add frame column
      if(vars_f[[2]] == "framest"){
        tmp_1 <- tmp_1 %>%
          dplyr::rename(frame = framest)
      }else{
        tmp_1 <- tmp_1 %>%
          dplyr::rename(frame = framesp)
      }

    }else{
      tmp_1 <- tmp %>%
        dplyr::group_by(sample,group,!!vars_f) %>%
        dplyr::summarise(avergae_exp = sum(norm_exp)/trans_numbers) %>%
        dplyr::filter(!!vars_f >= relative_distance[1] & !!vars_f <= relative_distance[2])
    }

    # check mode
    if(mode == "codon"){
      if(frame == TRUE){
        tmp_codon_positive <- tmp_1 %>%
          dplyr::filter(!!vars_f[[1]] >= 0) %>%
          dplyr::mutate(codon_pos = dplyr::if_else((!!vars_f[[1]] + 1)%% 3 == 0,(!!vars_f[[1]] + 1)/ 3,(!!vars_f[[1]] + 1)%/%3 + 1))

        tmp_codon_negtive <- tmp_1 %>%
          dplyr::filter(!!vars_f[[1]] < 0) %>%
          dplyr::mutate(codon_pos = dplyr::if_else(!!vars_f[[1]]%% 3 == 0,!!vars_f[[1]]/3,!!vars_f[[1]]%/%3 - 1))

        cb_df <- rbind(tmp_codon_positive,tmp_codon_negtive) %>%
          dplyr::group_by(sample,group,codon_pos,frame) %>%
          dplyr::summarise(avergae_exp = sum(avergae_exp)) %>%
          dplyr::mutate(pos = codon_pos)
      }else{
        tmp_codon_positive <- tmp_1 %>%
          dplyr::filter(!!vars_f >= 0) %>%
          dplyr::mutate(codon_pos = dplyr::if_else((!!vars_f + 1)%% 3 == 0,(!!vars_f + 1)/ 3,(!!vars_f + 1)%/%3 + 1))

        tmp_codon_negtive <- tmp_1 %>%
          dplyr::filter(!!vars_f < 0) %>%
          dplyr::mutate(codon_pos = dplyr::if_else(!!vars_f%% 3 == 0,!!vars_f/3,!!vars_f%/%3 - 1))

        cb_df <- rbind(tmp_codon_positive,tmp_codon_negtive) %>%
          dplyr::group_by(sample,group,codon_pos) %>%
          dplyr::summarise(avergae_exp = sum(avergae_exp)) %>%
          dplyr::mutate(pos = codon_pos)
      }
    }else{
      if(frame == TRUE){
        cb_df <- tmp_1 %>%
          dplyr::mutate(pos = !!vars_f[[1]])
      }else{
        cb_df <- tmp_1 %>%
          dplyr::mutate(pos = !!vars_f)
      }
    }

    # normalize
    total_sum_exp <- sum(cb_df$avergae_exp)
    ave_exp <- total_sum_exp/nrow(cb_df)
    cb_df$mean_exp <- cb_df$avergae_exp/ave_exp

    return(cb_df)
  }) -> df_final

  # ====================================================================================
  # plot
  # ====================================================================================
  gp <- unique(df_final$group)

  if(mode == "nt"){
    linewidth <- 0.1
  }else{
    linewidth <- 0.75
  }

  # deal with replicates
  condition <- !unique(is.na(gp))
  if(condition){
    df_plot <- df_final %>%
      dplyr::group_by(group,pos,mean_exp) %>%
      dplyr::summarise(mean_exp = mean(mean_exp)) %>%
      dplyr::mutate(group_name = "Metagene analysis")

    # line_layer <- geom_line(aes(x = pos,y = mean_exp))
    if(collapse == TRUE){
      line_layer <- do.call(geom_line,modifyList(list(mapping = aes(x = pos,y = mean_exp,color = group),
                                                      linewidth = linewidth),
                                                 geom_line_params))

      facet_layer <- do.call(facet_wrap,modifyList(list(facets = vars(group_name)),
                                                   facet_wrap_params))
    }else{
      line_layer <- do.call(geom_line,modifyList(list(mapping = aes(x = pos,y = mean_exp),
                                                      linewidth = linewidth),
                                                 geom_line_params))

      facet_layer <- do.call(facet_wrap,modifyList(list(facets = vars(group)),
                                                   facet_wrap_params))
    }

  }else{
    df_plot <- df_final %>%
      dplyr::mutate(group_name = "Metagene analysis")

    # check whether collapse plot
    if(collapse == TRUE){
      facet_layer <- do.call(facet_wrap,modifyList(list(facets = vars(group_name)),
                                                   facet_wrap_params))

      # line_layer <- geom_line(aes(x = pos,y = mean_exp,color = sample))
      line_layer <- do.call(geom_line,modifyList(list(mapping = aes(x = pos,y = mean_exp,color = sample),
                                                      linewidth = linewidth),
                                                 geom_line_params))
    }else{
      facet_layer <- do.call(facet_wrap,modifyList(list(facets = vars(sample)),
                                                   facet_wrap_params))

      # line_layer <- geom_line(aes(x = pos,y = mean_exp))
      if(frame == TRUE){
        line_layer <- do.call(geom_line,modifyList(list(mapping = aes(x = pos,y = mean_exp,
                                                                      color = factor(frame)),
                                                        linewidth = linewidth),
                                                   geom_line_params))
      }else{
        line_layer <- do.call(geom_line,modifyList(list(mapping = aes(x = pos,y = mean_exp),
                                                        linewidth = linewidth),
                                                   geom_line_params))
      }
    }

  }

  # plot
  ggplot(df_plot) +
    line_layer +
    theme_bw() +
    facet_layer +
    theme(panel.grid = element_blank(),
          strip.text = element_text(colour = "black",face = "bold",size = rel(1)),
          axis.text = element_text(colour = "black")) +
    xlab(paste("Distance from start/stop codon ","(",mode,")",sep = "")) +
    ylab("Average normalized footprint density\n (AU)")
}
