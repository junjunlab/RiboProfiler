globalVariables(c("abs_pos", "cds", "count", "pos", "sum_density" ,
                  "sum_pi", "tid", "utr5", "wi", "x", "y","rpm"))



#' Calculate Polarity Scores
#'
#' This function calculates polarity scores for gene expression data, filtered by CDS length.
#'
#' @param object ribosomeObj object.
#' @param minCounts Integer. Minimum number of counts required in the CDS region to include a gene. Default is 64.
#' @param norm_type The nomalization methods for ribosome density. "average" is calculated by
#' the count at each position divided by mean density across cds region. "rpm"
#' is calculated by the count at each position divided by the total counts and multiplied with 10^6.
#' Default is "average".
#' @param merge_rep Whether merge replicates. Default is FALSE.
#' @param group Character. A grouping variable for the output. Default is NULL, in which case the group is set to the current sample.
#' @param ... Useless args.
#'
#' @return A data frame with calculated polarity scores for each gene.
#'
#' @details
#' The function processes gene expression data by calculating a polarity score for each gene. It filters the data
#' by CDS length and counts, then computes the score based on the specified parameters.
#'
#' @examples
#' \dontrun{
#' calculatePolarity2(object = ribosomeObj)
#' }
#'
#' @export
setGeneric("calculate_polarity",
           function(object,
                    minCounts = 64,
                    norm_type = c("average","rpm"),
                    merge_rep = FALSE,
                    group = NULL, ...) standardGeneric("calculate_polarity"))





#' method for calculate_polarity
#'
#' @rdname calculate_polarity
#' @exportMethod calculate_polarity
setMethod("calculate_polarity",
          signature(object = "ribosomeObj"),
          function(object,
                   minCounts = 64,
                   norm_type = c("average","rpm"),
                   merge_rep = FALSE,
                   group = NULL,...){
            norm_type <- match.arg(norm_type,c("average","rpm"))

            # ==================================================================
            normed_file <- object@normalized.data

            # loop calculate polarity score
            sp <- unique(normed_file$sample)

            # x = 2
            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(normed_file,sample == sp[x])

              # gene annotation
              gene_anao <- object@longest.annotation
              colnames(gene_anao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                                       "cds_rg","exon_rg","utr5","cds","utr3")

              # filter cds length
              gene_lenth <- gene_anao %>%
                dplyr::select(gene_name,trans_id,utr5,cds,utr3)


              # filter counts in CDS region
              total_density <- tmp %>%
                dplyr::left_join(y = gene_lenth,by = "trans_id") %>%
                dplyr::filter(trans_pos >= utr5 & trans_pos <= utr5 + cds) %>%
                dplyr::group_by(trans_id) %>%
                dplyr::summarise(total_counts = sum(counts)) %>%
                dplyr::filter(total_counts >= minCounts,
                              sum_density = sum(!!rlang::ensym(norm_type)))

              # ========================================================================================
              # polarity calculation
              # ========================================================================================
              dt <- tmp %>%
                dplyr::filter(trans_id %in% total_density$trans_id) %>%
                dplyr::group_by(trans_id,trans_pos) %>%
                dplyr::summarise(density = sum(!!rlang::sym(norm_type))) %>%
                dplyr::left_join(y = gene_lenth,by = "trans_id") %>%
                dplyr::rename(pos = trans_pos)

              # calculate wi and pi for each gene
              if(is.null(group)){
                group <- sp[x]
              }

              ps_df <- dt |>
                dplyr::left_join(y = total_density,by = "trans_id",multiple = "all") |>
                dplyr::mutate(abs_pos = pos - utr5,
                              wi = (2*abs_pos -(cds + 1))/(cds - 1),
                              pi = density*wi/sum_density) |>
                dplyr::group_by(gene_name) |>
                dplyr::summarise(sum_pi = sum(pi)) |>
                dplyr::mutate(group = group,sample = sp[x],rep = tmp$rep[1])

              return(ps_df)
            }) -> ps

            # deal with replacates
            if(merge_rep == TRUE){
              ps <- ps |>
                dplyr::group_by(group,rep,gene_name) |>
                dplyr::summarise(sum_pi = mean(sum_pi)) |>
                dplyr::rename(sample = rep)
            }

            return(ps)
          }
)






#' Plot Polarity Scores
#'
#' This function creates a plot of polarity scores for genes based on RNA-seq data.
#' The function takes a data frame containing polarity scores and optionally a
#' grouping factor. It then creates a density plot of the polarity scores, with
#' the option to facet the plot by group.
#'
#' @param polarity_score_df A data frame containing the polarity scores.
#' @param aes_col_var A character string specifying the column variable for the
#'   aesthetic mapping. It can be either 'sample' or 'group'. Default is 'sample'.
#' @param linewidth A numeric value specifying the width of the density lines.
#'   Default is 1.
#' @param orf_col A character string specifying the color of the ORF line.
#'   Default is "grey".
#' @param orf_pos A numeric value specifying the position of the ORF line.
#'   Default is -0.2.
#' @param orf_width A numeric value specifying the width of the ORF line.
#'   Default is 4.
#' @param text_pos A numeric value specifying the position of the text label.
#'   Default is 0.1.
#' @param text_size A numeric value specifying the size of the text label.
#'   Default is 4.
#' @param facet A logical value specifying whether to facet the plot by group.
#'   Default is FALSE.
#' @param scales A character string specifying the scales for faceting.
#'   Default is "free".
#' @param ncol A numeric value specifying the number of columns for faceting.
#'   Default is 2.
#' @return A ggplot2 object representing the polarity score plot.
#' @import ggplot2
#' @export
polarity_score_plot <- function(polarity_score_df = NULL,
                                aes_col_var = c("sample","group"),
                                linewidth = 1,
                                orf_col = "grey",
                                orf_pos = -0.2,
                                orf_width = 4,
                                text_pos = 0.1,
                                text_size = 4,
                                facet = FALSE,
                                scales = "free",
                                ncol = 2){
  # check args
  aes_col_var <- match.arg(aes_col_var,c("sample","group"))

  # line layer
  line_layer <- geom_density(aes(color = .data[[aes_col_var]]),linewidth = linewidth)

  # plot
  pmain <-
    ggplot(polarity_score_df,aes(x = sum_pi,y = after_stat(density))) +
    line_layer +
    # geom_density() +
    theme_bw() +
    theme(axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold.italic",size = rel(1))) +
    scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1),
                       limits = c(-1,1),
                       labels = c("-1","-0.5","0","0.5","1")) +
    ylab("Relative gene density") +
    xlab("Polarity score") +
    # annotation
    geom_line(data = data.frame(x = c(-1,1),y = rep(orf_pos,2)),
              aes(x = x,y = y),
              color = orf_col,
              linewidth = orf_width) +
    geom_text(data = data.frame(x = 0,y = text_pos),
              aes(x = x,y = y,label = "ORF (5'->3')"),
              fontface = "bold",size = text_size)

  # check whether facet for plot
  if(facet == TRUE){
    pmain +
      facet_wrap(~group,scales = scales,ncol = ncol)
  }else{
    pmain
  }
}
