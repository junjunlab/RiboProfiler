#' Create a footprint heatmap
#'
#' This function creates a footprint heatmap for ribo-seq data. The heatmap
#' visualizes the count of reads at different distances from the start
#' or stop codon relative to the 5' end of the transcript. The type of
#' distance can be specified as 'relst' for the distance from the start
#' codon or 'relsp' for the distance from the stop codon.
#'
#' @param object ribosomeObj object.
#' @param type A character string specifying the type of distance to use.
#'   It can be either 'relst' or 'relsp'. Default is 'relst'.
#' @param rel_distance A numeric vector specifying the relative distance
#'   from the 5' end to the start or stop codon. Default is NULL.
#' @param read_length A numeric vector specifying the read length range.
#'   Default is c(20,35).
#' @param tile_border A character string specifying the color of the tile
#'   border. Default is "grey".
#' @param bar_col A character string specifying the color of the bar.
#'   Default is "grey".
#' @param bar_width A numeric value specifying the width of the bar.
#' @param scale_x X panel size. Default is 0.1.
#' @param scale_y Y panel size. Default is 0.1.
#' @param xline Whether show x side plot with line plot. Default is FASLE.
#' @param xline_col The color of x side line when xline = TRUE. Default is "black".
#'   Default is 0.75.
#' @param log Whether do log transformation for counts. Default is FALSE.
#' @param xseq The x axis interval. Default is 3.
#' @param ... Useless args.
#'
#' @return A ggplot2 object representing the footprint heatmap.
#' @importFrom dplyr filter group_by summarise
#' @import ggplot2
#' @importFrom ggside geom_xsidecol scale_xsidey_continuous geom_ysidecol
#' @export
setGeneric("footprint_heatmap",
           function(object,
                    type = c("relst","relsp"),
                    rel_distance = NULL,
                    read_length = c(20,35),
                    tile_border = "grey",
                    bar_col = "grey",
                    bar_width = 0.75,
                    scale_x = 0.1,
                    scale_y = 0.1,
                    xseq = 3,
                    xline = FALSE,
                    xline_col = "black",
                    log = FALSE, ...) standardGeneric("footprint_heatmap"))







#' method for footprint_heatmap
#'
#' @rdname footprint_heatmap
#' @exportMethod footprint_heatmap
setMethod("footprint_heatmap",
          signature(object = "ribosomeObj"),
          function(object,
                   type = c("relst","relsp"),
                   rel_distance = NULL,
                   read_length = c(20,35),
                   tile_border = "grey",
                   bar_col = "grey",
                   bar_width = 0.75,
                   scale_x = 0.1,
                   scale_y = 0.1,
                   xseq = 3,
                   xline = FALSE,
                   xline_col = "black",
                   log = FALSE,...){
            # ==========================================================================
            # check args
            # ==========================================================================
            type <- match.arg(type,c("relst","relsp"))

            if(type == "relst"){
              if(is.null(rel_distance)){
                rel_pos <- c(-18,12)
              }else{
                rel_pos <- rel_distance
              }

              xlab <- "Distance of 5'end to start codon (nt)"
            }else{
              if(is.null(rel_distance)){
                rel_pos <- c(-25,5)
              }else{
                rel_pos <- rel_distance
              }

              xlab <- "Distance of 5'end to stop codon (nt)"
            }

            # filter read length
            qc_data <- object@raw.counts
            qc_data <- qc_data |>
              dplyr::mutate(length = as.numeric(as.character(len))) |>
              dplyr::filter(length >= read_length[1] & length <= read_length[2])

            # ==========================================================================
            # process data
            # ==========================================================================
            ht <- qc_data |>
              dplyr::filter(.data[[type]] >= rel_pos[1] & .data[[type]] <= rel_pos[2]) |>
              dplyr::group_by(sample,length,.data[[type]]) |>
              dplyr::summarise(count = sum(counts))

            if(log == TRUE){
              ht$count <- log2(ht$count + 1)
            }


            ht_x_df <- ht |>
              dplyr::group_by(sample,.data[[type]]) |>
              dplyr::summarise(count = sum(count))

            ht_y_df <- ht |>
              dplyr::group_by(sample,length) |>
              dplyr::summarise(count = sum(count))

            # ==========================================================================
            # plot
            # ==========================================================================
            # check xlayer
            if(xline == TRUE){
              xlayer <- ggside::geom_xsideline(data = ht_x_df,
                                               aes(x = .data[[type]],y = count),color = xline_col)
            }else{
              xlayer <- ggside::geom_xsidecol(data = ht_x_df,
                                              aes(x = .data[[type]],y = count),fill = bar_col,
                                              width = bar_width)
            }

            # plot
            ggplot(ht) +
              geom_tile(aes(x = .data[[type]],y = length,fill = count),color = tile_border) +
              # x barplot
              xlayer +
              ggside::scale_xsidey_continuous(breaks = NULL) +
              # y barplot
              ggside::geom_ysidecol(data = ht_y_df,
                                    aes(x = count,y = length),fill = bar_col,
                                    width = bar_width) +
              ggside::scale_ysidex_continuous(breaks = NULL) +
              # others
              theme_bw() +
              scale_fill_viridis_c(option = "mako",direction = -1) +
              scale_y_continuous(breaks = seq(read_length[1],read_length[2],2),
                                 labels = seq(read_length[1],read_length[2],2)) +
              scale_x_continuous(breaks = seq(rel_pos[1],rel_pos[2],xseq),
                                 labels = seq(rel_pos[1],rel_pos[2],xseq)) +
              coord_cartesian(expand = 0) +
              theme(axis.text = element_text(colour = "black"),
                    panel.grid = element_blank(),
                    strip.text = element_text(face = "bold.italic",size = rel(1)),
                    ggside.panel.background = element_blank(),
                    ggside.panel.border = element_blank(),
                    ggside.panel.scale.x = scale_x,
                    ggside.panel.scale.y = scale_y) +
              xlab(xlab) +
              ylab("Read length (nt)") +
              facet_wrap(~sample,scales = "free")

          }
)
