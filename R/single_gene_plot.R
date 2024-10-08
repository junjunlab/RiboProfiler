globalVariables(c("RPM"))

#' Single Gene Plot
#'
#' This method generates a plot showing ribosome positioning along a transcript for a specific gene.
#'
#' @param object An object of class `ribosomeObj` containing normalized ribosome
#' profiling data and gene annotations.
#' @param select_gene A character string specifying the gene of interest to plot.
#' Must correspond to a gene name in the provided object.
#' @param merge_rep Logical. If `TRUE`, replicates will be merged by taking the mean RPM. Default is `FALSE`.
#' @param mode Character. Determines the unit of x-axis as "nt" for nucleotides
#' or "codon" for codons/amino acids. Default is `"nt"`.
#' @param structure_position Character. Defines the position of transcript
#' structure ("top" or "bottom"). Default is `"top"`.
#' @param orf_col Color code for the coding sequence (CDS) region in the plot. Default is `"grey65"`.
#' @param utr_col Color code for the untranslated regions (UTRs) in the plot. Default is `"grey"`.
#' @param utr_width Numeric. The line width for UTRs in the plot. Default is `1`.
#' @param cds_width Numeric. The line width for CDS in the plot. Default is `8`.
#' @param cds_region_scale_x Numeric. Scales the x-axis for the CDS region. Default is `0.1`.
#' @param cds_region_anno Data frame. An optional annotation for the CDS region to be plotted.
#' @param cds_anno_col Color code for CDS annotations. Default is `"#996600"`.
#' @param anno_region_width Numeric. The width for the annotated regions. Default is `5`.
#' @param geom_col_params List. Additional parameters to be passed to `geom_col`
#' for customizing the bar plot aesthetics.
#' @param ... Other arguments passed to internal plotting functions.
#'
#' @return A ggplot object displaying the ribosome footprint density across the selected gene transcript.
#'
#' @examples
#' \dontrun{# Assume 'ribosome_data' is a preloaded object of class 'ribosomeObj'
#' single_gene_plot(ribosome_data, select_gene = "GCN4", mode = "codon", merge_rep = TRUE)}
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom ggside geom_xsidesegment
#' @export
setGeneric("single_gene_plot",
           function(object,
                    select_gene = NULL,
                    merge_rep = FALSE,
                    mode = c("nt","codon"),
                    structure_position = c("top","bottom"),
                    orf_col = "grey65",
                    utr_col = "grey",
                    utr_width = 1,
                    cds_width = 8,
                    cds_region_scale_x = 0.1,
                    cds_region_anno = NULL,
                    cds_anno_col = "#996600",
                    anno_region_width = 5,
                    geom_col_params = list(),
                    ...) standardGeneric("single_gene_plot"))






#' method for single_gene_plot
#'
#' @rdname single_gene_plot
#' @exportMethod single_gene_plot
setMethod("single_gene_plot",
          signature(object = "ribosomeObj"),
          function(object,
                   select_gene = NULL,
                   merge_rep = FALSE,
                   mode = c("nt","codon"),
                   structure_position = c("top","bottom"),
                   orf_col = "grey65",
                   utr_col = "grey",
                   utr_width = 1,
                   cds_width = 8,
                   cds_region_scale_x = 0.1,
                   cds_region_anno = NULL,
                   cds_anno_col = "#996600",
                   anno_region_width = 5,
                   geom_col_params = list(),
                   ...){
            mode <- match.arg(mode,c("nt","codon"))
            structure_position <- match.arg(structure_position,c("top","bottom"))
            # ==================================================================
            # prepare data
            # ==================================================================
            norm_df <- object@normalized.data

            # "THO1" "GCN4"
            ref <- object@longest.annotation %>%
              dplyr::filter(gene_name == select_gene)

            xlimit <- c(-(ref$utr5 - 1),ref$cds + ref$utr3)
            xlabel <- "Ribosome along transcript position\n (nucleotides / nts)"

            # subset gene
            g_df <- subset(norm_df,trans_id == ref$tid) %>%
              dplyr::group_by(sample,rep,group,trans_id,trans_pos) %>%
              dplyr::summarise(RPM = sum(rpm),.groups = "drop") %>%
              # dplyr::filter(trans_pos > 0) %>%
              dplyr::left_join(y = ref,by = c("trans_id" = "tid")) %>%
              dplyr::mutate(pos = trans_pos - utr5)

            # check mode
            if(mode == "codon"){
              g_df <- g_df %>%
                dplyr::filter(pos > 0) %>%
                dplyr::mutate(pos = ifelse(pos %% 3 == 0,pos/3,pos%/%3 + 1)) %>%
                dplyr::group_by(sample,rep,group,gene_name,trans_id,pos) %>%
                dplyr::summarise(RPM = sum(RPM),.groups = "drop")

              xleft <- ifelse(-(ref$utr5 - 1) %% 3 == 0,-(ref$utr5 - 1)/3,-(ref$utr5 - 1)%/%3 + 1)
              xright <- ifelse((ref$cds + ref$utr3) %% 3 == 0,(ref$cds + ref$utr3)/3,(ref$cds + ref$utr3)%/%3 + 1)
              xlimit <- c(xleft,xright)

              xlabel <- "Ribosome along transcript position\n (codons / amino acids)"
            }

            # ==================================================================
            # deal with replicates
            # ==================================================================
            if(merge_rep == TRUE){
              g_df <- g_df %>%
                dplyr::group_by(rep,group,gene_name,trans_id,pos) %>%
                dplyr::summarise(RPM = mean(RPM),.groups = "drop")
            }

            # ==================================================================
            # main plot
            # ==================================================================
            pmian <-
              ggplot(g_df) +
              # geom_col(aes(x = pos,y = RPM),width = 1) +
              do.call(geom_col,modifyList(list(mapping = aes(x = pos,y = RPM),
                                               width = 1),
                                          geom_col_params)) +
              facet_grid(cols = vars(gene_name),rows = vars(sample)) +
              theme_bw() +
              theme(axis.text = element_text(colour = "black"),
                    panel.grid = element_blank(),
                    strip.placement = "outside",
                    strip.text = element_text(face = "bold",size = rel(1)),
                    strip.text.x = element_text(face = "bold.italic",size = rel(1)),
                    strip.background = element_blank(),
                    strip.text.y.right = element_text(angle = 0,hjust = 0)) +
              xlim(xlimit) +
              xlab(xlabel) +
              ylab("Ribosome footprint density (RPM)")

            # ==================================================================
            # gene structure
            # ==================================================================
            trans_struc <- data.frame(x = xlimit[1],xend = xlimit[2],y = 1,yend = 1)
            if(mode == "codon"){
              cds_struc <- data.frame(x = 1,
                                      xend = ifelse(ref$cds%%3==0,ref$cds/3,ref$cds%/%3 + 1),
                                      y = 1,yend = 1)
            }else{
              cds_struc <- data.frame(x = 1,xend = ref$cds,y = 1,yend = 1)
            }

            # whether add new structure
            if(!is.null(cds_region_anno)){
              tmp_struc <- subset(cds_region_anno,gene_name == select_gene)
              struc_anno <- data.frame(x = as.numeric(tmp_struc$start),
                                       xend = as.numeric(tmp_struc$end),
                                       y = 1,yend = 1)


              struc_anno_layer <- ggside::geom_xsidesegment(data = struc_anno,
                                                            mapping = aes(x = x,xend = xend,y = y,yend = yend),
                                                            linewidth = anno_region_width,
                                                            color = cds_anno_col)
            }else{
              struc_anno_layer <- NULL
            }

            # plot
            pmian +
              ggside::geom_xsidesegment(data = trans_struc,
                                        mapping = aes(x = x,xend = xend,y = y,yend = yend),
                                        linewidth = utr_width,color = utr_col) +
              ggside::geom_xsidesegment(data = cds_struc,
                                        mapping = aes(x = x,xend = xend,y = y,yend = yend),
                                        linewidth = cds_width,color = orf_col) +
              struc_anno_layer +
              theme(ggside.panel.background = element_blank(),
                    ggside.panel.border = element_blank(),
                    ggside.panel.scale.x = cds_region_scale_x) +
              ggside::scale_xsidey_continuous(breaks = NULL) +
              ggside::ggside(x.pos = structure_position,collapse = "x")

          })
