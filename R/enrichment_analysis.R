globalVariables(c("compare_group"))

#' Compute Sliding Window Average
#'
#' This function computes the sliding window average of a numeric vector. The window averages
#' are calculated for each point where the center of the window is moved stepwise
#' across the vector. The window has an adjustable size, and the step size can also be customized.
#'
#' @param ratio A numeric vector containing the data for which the sliding window average is to be calculated.
#' @param window An integer specifying the size of the sliding window. Default is \code{7}.
#' @param step An integer specifying the step size for moving the window across the vector. Default is \code{1}.
#'
#' @details
#' The sliding window average is computed by taking the mean of values within a specified
#' window size around each central value. For the boundary regions (near the start and end), the function
#' retains the original values without applying the moving average.
#'
#' The cumulative sum of the input vector (\code{ratio}) is computed for efficient window sum calculation.
#' Boundary handling ensures that the first and last `window-1` data points remain unaltered.
#'
#' @return
#' Returns a numeric vector of the same length as the input \code{ratio}, where each value corresponds
#' to the sliding window average. Boundary regions (near the start and end) retain their original values.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' result <- slide_window_average(ratio = data, window = 3, step = 1)
#' print(result)
#' }
#'
#' @export
slide_window_average <- function(ratio = NULL, window = 7, step = 1) {
  # start <- window - 1
  start <- (window - 1)/2 + 1
  winLen <- length(ratio)
  tmp_data <- numeric(winLen)

  # exclude out boundary
  tmp_data[1:start] <- ratio[1:start]
  tmp_data[(winLen - start + 1):winLen] <- ratio[(winLen - start + 1):winLen]

  # cusums
  cum_sum <- c(0, cumsum(ratio))

  # moving average
  for (j in seq(start, winLen - start, by = step)) {
    left <- max(1, j - floor((window - 1) / 2))
    right <- min(winLen, j + floor((window - 1) / 2))
    window_sum <- cum_sum[right + 1] - cum_sum[left]
    tmp_data[j] <- window_sum / (right - left + 1)
  }

  return(tmp_data)
}







#' Plot Sequence Frame
#'
#' This function plots the sequence frame of a given CDS (Coding DNA Sequence) with
#' annotations for start and stop codons. The function can be used to visualize the
#' positions of these codons relative to the reading frame of the transcript.
#'
#' @param cds_seq A character string or Biostrings object representing the coding DNA sequence (CDS).
#' @param seq_mark A character string specifying the type of codon to annotate:
#'   \itemize{
#'     \item \code{"both"} for both start and stop codons.
#'     \item \code{"start"} for start codons only.
#'     \item \code{"stop"} for stop codons only.
#'   }
#'   Default is \code{"both"}.
#' @param trans_id A character string or integer specifying the transcript ID for which the CDS sequence should be plotted.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{plot}: A ggplot2 object representing the sequence frame plot.
#'     \item \code{df_pos}: A data frame containing the positions of start and stop codons in the sequence.
#'   }
#'
#'
#' @export
seq_frame_plot <- function(cds_seq = NULL,
                           seq_mark = c("both","start","stop"),
                           trans_id = NULL){
  seq_mark <- match.arg(seq_mark,c("both","start","stop"))

  # load sequence
  cds_fa <- Biostrings::readDNAStringSet(cds_seq)
  names(cds_fa) <- sapply(strsplit(names(cds_fa),split = " "),"[",1)
  seq <- cds_fa[trans_id]

  # find ATG/TAG/TGA/TAA positions
  aug_seq_pos <- stringr::str_locate_all(as.character(seq),pattern = "ATG") %>%
    data.frame() %>%
    dplyr::mutate(frame = (start - 1)%%3,type = "start_codon")

  uxx_seq_pos <- stringr::str_locate_all(as.character(seq),pattern = "TAG|TGA|TAA") %>%
    data.frame() %>%
    dplyr::mutate(frame = (start - 1)%%3,type = "stop_codon")

  # check mark
  if(seq_mark == "start"){
    df_pos <- aug_seq_pos
  }else if(seq_mark == "stop"){
    df_pos <- uxx_seq_pos
  }else{
    df_pos <- rbind(aug_seq_pos,uxx_seq_pos)
  }

  # plot
  p <-
    ggplot(df_pos) +
    geom_rect(aes(xmin = start,xmax = end,ymin = 0,ymax = 1,
                  fill = type),show.legend = F) +
    facet_wrap(~frame,ncol = 1,strip.position = "left") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.text.y.left = element_text(angle = 0),
          panel.spacing.y = unit(1,"mm"),
          plot.margin = margin(t = 0,unit = "cm")) +
    scale_fill_manual(values = c(start_codon = "#CC3333",stop_codon = "black"),
                      name = "")

  return(list(plot = p,df_pos = df_pos))
}










# ==============================================================================
# methods
# ==============================================================================

#' Enrichment Analysis for Ribosome Profiling Data
#'
#' This function performs enrichment analysis on ribosome profiling data contained in
#' a \code{ribosomeObj} object. The analysis compares total ribosome counts with
#' seRP (site-specific ribosome profiling) counts across various metrics and aggregates
#' the results based on codon or nucleotide positions.
#'
#' @param object An object of class \code{ribosomeObj} which contains normalized data and gene annotation.
#' @param total A character vector specifying the sample names for total ribosome.
#' @param seRP A character vector specifying the sample names for seRP, corresponding to `total`.
#' @param mode A character string specifying the analysis mode. Possible values are:
#'   \itemize{
#'     \item \code{"codon"}: Perform analysis at the codon level.
#'     \item \code{"nt"}: Perform analysis at the nucleotide level.
#'   }
#'   Default is \code{"codon"}.
#' @param rep_name A character vector specifying replicate names for each sample.
#' @param compare_group_name A character vector specifying the group names for comparison.
#' @param ... Additional arguments (not currently used).
#'
#' @details
#' This function performs a sliding window analysis of ribosome profiling data, where codon
#' occupancy or nucleotide occupancy is calculated based on combined seRP and total ribosome counts.
#' The results are summarized in terms of the ratio \code{seRP / total}.
#' The output is stored in the \code{enriched.data} slot of the \code{ribosomeObj} object.
#'
#' The function processes the data by:
#' \enumerate{
#'   \item Merging normalized data with gene annotations by transcript ID.
#'   \item Converting transcript positions based on the selected mode (codon or nucleotide).
#'   \item Grouping the data by \code{sample}, \code{gene_name}, \code{trans_id}, \code{type}, and position.
#'   \item Calculating the seRP/total ratio for each group.
#'   \item Aggregating the results to generate mean/seRP ratios by gene and position.
#' }
#'
#' The final result is added to the \code{enriched.data} slot of the \code{ribosomeObj} object,
#' and the unit attribute (codon or nucleotide) is attached to the results.
#'
#' @return A modified \code{ribosomeObj} object with the enrichment results added to
#' the \code{enriched.data} slot.
#'
#' @examples
#' # Assuming riboObj is a ribosomeObj object with normalized data and annotations:
#' \dontrun{
#' riboObj <- enrichment_analysis(
#'   object = riboObj,
#'   total = c("sample_total_1", "sample_total_2"),
#'   seRP = c("sample_seRP_1", "sample_seRP_2"),
#'   mode = "codon",
#'   rep_name = c("replicate_1", "replicate_2"),
#'   compare_group_name = c("groupA", "groupB")
#' )}
#'
#' @export
setGeneric("enrichment_analysis",
           function(object,
                    total = NULL,
                    seRP = NULL,
                    mode = c("codon","nt"),
                    rep_name = NULL,
                    compare_group_name = NULL, ...) standardGeneric("enrichment_analysis"))







#' method for enrichment_analysis
#'
#' @rdname enrichment_analysis
#' @exportMethod enrichment_analysis
setMethod("enrichment_analysis",
          signature(object = "ribosomeObj"),
          function(object = NULL,
                   total = NULL,
                   seRP = NULL,
                   mode = c("codon","nt"),
                   rep_name = NULL,
                   compare_group_name = NULL,...){
            mode <- match.arg(mode,c("codon","nt"))
            # =====================================================================
            normed_df_anno <- object@normalized.data %>%
              dplyr::mutate(type = dplyr::if_else(sample %in% total,"total","seRP"))

            # load geneinfo
            geneInfo <- object@longest.annotation
            colnames(geneInfo) <- c('id','gene_name','gene_id','trans_id','chr','strand',
                                    'cds_region','exon_region','utr5','cds','utr3')

            geneInfo <- geneInfo %>%
              dplyr::select(gene_name,trans_id,utr5,cds,utr3)

            # x = 1
            purrr::map_df(seq_along(total),function(x){
              compare_var <- c(total[x],seRP[x])

              tmp <- normed_df_anno %>%
                dplyr::filter(sample %in% compare_var) %>%
                dplyr::left_join(y = geneInfo,by = "trans_id") %>%
                dplyr::mutate(trans_pos = trans_pos - utr5) %>%
                dplyr::filter(trans_pos > 0 & trans_pos <= cds)


              if(mode == "codon"){
                tmp2 <- tmp %>%
                  dplyr::mutate(pos = dplyr::if_else(trans_pos%% 3 == 0,trans_pos/3,trans_pos%/%3 + 1)) %>%
                  dplyr::group_by(sample,gene_name,trans_id,type,pos) %>%
                  dplyr::summarise(exp = sum(norm_exp),.groups = "drop")
              }else{
                tmp2 <- tmp %>%
                  dplyr::mutate(pos = trans_pos) %>%
                  dplyr::group_by(sample,gene_name,trans_id,type,pos) %>%
                  dplyr::summarise(exp = sum(norm_exp),.groups = "drop")
              }

              # calculate ratio of seRP/total
              tmp3 <- tmp2 %>%
                dplyr::select(-sample) %>%
                tidyr::spread(type,exp) %>%
                dplyr::mutate(ratio = .data[["seRP"]]/.data[["total"]]) %>%
                tidyr::replace_na(list(ratio = 0)) %>%
                dplyr::mutate(compare_group = compare_group_name[x],
                              rep = rep_name[x])

              return(tmp3)
            }) -> enrich_df

            # calculate sd
            data.table::setDT(enrich_df)
            enrich_df <- enrich_df[, .(ratio = mean(ratio)),
                                   by = .(gene_name, trans_id, pos, compare_group,rep)]

            # assign a unit attribute
            attr(enrich_df,"unit") <- mode

            object@enriched.data <- enrich_df

            return(object)
          }
)







#' Generate Ribosome Enrichment Track Plot
#'
#' This function creates enrichment track plots to visualize ribosome profiling data.
#' Users can visualize either ribosome footprint density or the enrichment ratio
#' (IP/Total) across transcripts. The function generates highly customizable plots,
#' including options for codon or nucleotide-level resolution, plotting as line or
#' column, and more.
#'
#' @param object An object of class \code{ribosomeObj}, the ribosome profiling object
#'   that contains gene annotations, normalized data, and enriched data.
#' @param select_gene A character vector specifying the gene name(s) to plot.
#' @param select_sample A character vector specifying the sample name(s) to plot.
#'   Default is \code{NULL}, which plots all available samples.
#' @param assay A character string specifying whether to plot ribosome footprint density
#'   or enrichment data. Must be one of \code{"normalized"} or \code{"enriched"}. Default is \code{"normalized"}.
#' @param geom_type A character string defining the geometry type for the plot. Can be
#'   either \code{"line"} (line plot) or \code{"col"} (column/bar plot). Default is \code{"line"}.
#' @param merge_rep Logical; if \code{TRUE}, data from replicates are merged into one plot. Default is \code{TRUE}.
#' @param collapse Logical; if \code{TRUE}, replicates or sample groups are collapsed into a single plot. Default is \code{FALSE}.
#' @param sample_color A named vector of colors for samples. Default is \code{NULL}, which uses ggplot's default color palette.
#' @param average_smooth Logical; if \code{TRUE}, a sliding window average is applied to the ratio data. Default is \code{TRUE}.
#' @param window_size Integer specifying the window size for the sliding average. Default is \code{7}.
#' @param mode A character string specifying the unit for plotting. Can be either \code{"codon"} or \code{"nt"} (nucleotide). Default is \code{"codon"}.
#' @param add_frame_anno Logical; if \code{TRUE}, frame annotations are added to the plot for nucleotide plots. Default is \code{FALSE}.
#' @param frame_relheight Numeric value specifying the relative height of the frame annotation plot. Default is \code{0.1}.
#' @param seq_mark A character string specifying the type of codon to annotate. Must be one of \code{"both"}, \code{"start"}, or \code{"stop"}. Default is \code{"both"}.
#' @param cds_region_anno Optional data frame containing coding sequence annotation.
#'   Default is \code{NULL}.
#' @param cds_col A character string specifying the color for visualizing CDS regions. Default is \code{"grey"}.
#' @param cds_anno_col A character string specifying the color for coding sequence annotation. Default is \code{"#996600"}.
#' @param cds_region_width Numeric value specifying the width of the CDS structure annotation. Default is \code{5}.
#' @param cds_region_scale_x Numeric value for scaling the coding sequence on the x-axis. Default is \code{NULL}.
#' @param geom_line_list A list of ggplot2 aesthetic modifications for the line geometry. Default is an empty list.
#' @param geom_col_list A list of ggplot2 aesthetic modifications for the column geometry. Default is an empty list.
#' @param facet_grid_list A list of parameters for modifying `ggplot2::facet_grid()` layout. Default is an empty list.
#' @param ... Further arguments to be passed to internal functions.
#'
#' @return Returns a \code{ggplot} object or a combined plot list if multiple genes are plotted.
#'
#' @details
#' The function visualizes ribosome profiling data by plotting it as either ribosome footprint
#' density (from \code{object@normalized.data}) or enrichment ratio (from \code{object@enriched.data}.
#'
#' Users can choose to display data either by codons or nucleotide, select start and stop codons for annotation,
#' apply optional smoothing using a sliding window average, and add customizable gene structure annotations.
#'
#' If `add_frame_anno = TRUE`, frame annotations will be added for a nucleotide plot.
#'
#' @examples
#' \dontrun{
#' data <- ribosomeObj # Assuming ribosomeObj is an object of class ribosomeObj
#' p <- enrichment_track_plot(
#'   object = data,
#'   select_gene = "GeneA",
#'   select_sample = c("Sample1", "Sample2"),
#'   assay = "normalized",
#'   geom_type = "line",
#'   merge_rep = TRUE
#' )
#' print(p) # Visualize the plot
#' }
#'
#' @export
setGeneric("enrichment_track_plot",
           function(object,
                    select_gene = NULL,
                    select_sample = NULL,
                    assay = c("normalized","enriched"),
                    geom_type = c("line","col"),
                    merge_rep = TRUE,
                    collapse = FALSE,
                    sample_color = NULL,
                    average_smooth = TRUE,
                    window_size = 7,
                    mode = c("codon","nt"),
                    add_frame_anno = FALSE,
                    frame_relheight = 0.1,
                    seq_mark = c("both","start","stop"),
                    cds_region_anno = NULL,
                    cds_col = "grey",
                    cds_anno_col = "#996600",
                    cds_region_width = 5,
                    cds_region_scale_x = NULL,
                    geom_line_list = list(),
                    geom_col_list = list(),
                    facet_grid_list = list(), ...) standardGeneric("enrichment_track_plot"))







#' method for enrichment_track_plot
#'
#' @rdname enrichment_track_plot
#' @exportMethod enrichment_track_plot
setMethod("enrichment_track_plot",
          signature(object = "ribosomeObj"),
          function(object,
                   select_gene = NULL,
                   select_sample = NULL,
                   assay = c("normalized","enriched"),
                   geom_type = c("line","col"),
                   merge_rep = TRUE,
                   collapse = FALSE,
                   sample_color = NULL,
                   average_smooth = TRUE,
                   window_size = 7,
                   mode = c("codon","nt"),
                   add_frame_anno = FALSE,
                   frame_relheight = 0.1,
                   seq_mark = c("both","start","stop"),
                   cds_region_anno = NULL,
                   cds_col = "grey",
                   cds_anno_col = "#996600",
                   cds_region_width = 5,
                   cds_region_scale_x = NULL,
                   geom_line_list = list(),
                   geom_col_list = list(),
                   facet_grid_list = list(),...){
            assay <- match.arg(assay,c("normalized","enriched"))
            mode <- match.arg(mode,c("codon","nt"))
            seq_mark <- match.arg(seq_mark,c("both","start","stop"))
            geom_type <- match.arg(geom_type,c("line","col"))
            # ============================================================================
            # gene annotation
            # ============================================================================
            ganao <- object@longest.annotation
            colnames(ganao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                                 "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")
            ganao <- ganao[,c("gene_name","gene_id","trans_id","5UTR_length","CDS_length")]
            ganao <- subset(ganao,gene_name %in% select_gene)

            # ============================================================================
            # load enriched data
            # ============================================================================
            if(assay == "enriched"){
              ylabel <- "Enrichment (IP / Total)"
              ratio_df <- object@enriched.data

              if(is.null(select_sample)){
                select_sample <- unique(ratio_df$compare_group)
              }

              ratio_df_new <- ratio_df %>%
                dplyr::filter(gene_name %in% select_gene & compare_group %in% select_sample)
            }else{
              ylabel <- "Ribosome footprint density"
              ratio_df <- object@normalized.data

              if(is.null(select_sample)){
                select_sample <- unique(ratio_df$sample)
              }

              # transform to codon pos
              ratio_df_new <- ratio_df %>%
                dplyr::select(trans_pos,trans_id,sample,rep,group,norm_exp) %>%
                dplyr::left_join(y = ganao,by = "trans_id") %>%
                dplyr::filter(gene_name %in% select_gene & sample %in% select_sample) %>%
                dplyr::mutate(cdsft = dplyr::if_else(CDS_length%%3 == 0,1,0)) %>%
                # filter cds length %%3 == 0
                dplyr::filter(cdsft == 1) %>%
                dplyr::mutate(trans_pos = trans_pos - `5UTR_length`)

              # check unit
              if(mode == "codon"){
                ratio_df_new <- ratio_df_new %>%
                  # trans_pos into codon pos
                  dplyr::mutate(pos = dplyr::if_else(trans_pos%%3 == 0, trans_pos/3,trans_pos%/%3 + 1)) %>%
                  dplyr::filter(pos > 0 & pos < CDS_length) %>%
                  dplyr::group_by(sample,rep,group,gene_name,trans_id,pos) %>%
                  dplyr::summarise(ratio = mean(norm_exp),.groups = "drop") %>%
                  dplyr::rename(compare_group = sample)
              }else{
                ratio_df_new <- ratio_df_new %>%
                  dplyr::filter(trans_pos > 0 & trans_pos < CDS_length) %>%
                  dplyr::rename(pos = trans_pos) %>%
                  dplyr::group_by(sample,rep,group,gene_name,trans_id,pos) %>%
                  dplyr::summarise(ratio = mean(norm_exp),.groups = "drop") %>%
                  dplyr::rename(compare_group = sample)
              }
            }

            # ============================================================================
            # generate track plot
            # ============================================================================

            # loop plot
            # x = 1
            lapply(seq_along(select_gene),function(x){
              if(x != length(select_gene)){
                strip_text_y <- element_blank()
              }else{
                strip_text_y <- element_text(angle = 0,hjust = 0)
              }

              gene_df <- subset(ganao,gene_name == select_gene[x])

              if(mode == "codon"){
                cds_len <- gene_df$CDS_length/3
              }else{
                cds_len <- gene_df$CDS_length
              }


              # ================================================================
              # line layer
              if(merge_rep == TRUE){
                col_var <- rlang::sym("rep")
              }else{
                col_var <- rlang::sym("compare_group")
              }

              # check geom type
              if(geom_type == "line"){
                geom_plot_layer <- do.call(geom_line,
                                           modifyList(list(mapping = aes(x = pos,y = ratio,
                                                                         color = !!col_var),
                                                           show.legend = FALSE),
                                                      geom_line_list))
              }else{
                geom_plot_layer <- do.call(geom_col,
                                           modifyList(list(mapping = aes(x = pos,y = ratio,
                                                                         fill = !!col_var),
                                                           show.legend = FALSE),
                                                      geom_col_list))
              }

              df_plotp <- ratio_df_new %>% dplyr::filter(gene_name == select_gene[x])

              # ================================================================
              # whether do average smooth
              # ================================================================
              if(average_smooth == TRUE){
                sp <- unique(df_plotp$compare_group)

                # x = 1
                purrr::map_df(seq_along(sp),function(x){
                  df_rg <- data.frame(pos = 1:cds_len)

                  tmp <- subset(df_plotp,compare_group == sp[x])

                  tmp1 <- df_rg %>%
                    dplyr::left_join(y = tmp[,c("pos","ratio")],by = "pos") %>%
                    tidyr::replace_na(list(ratio = 0))

                  smoothed_ratio <- slide_window_average(ratio = tmp1$ratio,window = window_size)

                  tmp_out <- data.frame(gene_name = tmp$gene_name[1],
                                        trans_id = tmp$trans_id[1],
                                        pos = tmp1$pos,
                                        compare_group = tmp$compare_group[1],
                                        rep = tmp$rep[1],
                                        ratio = smoothed_ratio)

                }) -> df_plotp_smoothed
              }else{
                df_plotp_smoothed <- df_plotp
              }

              # ================================================================
              # whether merge replicates
              # ================================================================
              if(merge_rep == TRUE){
                df_plotp_merged_rep <- df_plotp_smoothed %>%
                  dplyr::group_by(rep,gene_name,trans_id,pos) %>%
                  dplyr::summarise(ratio = mean(ratio),.groups = "drop")

                row_var <- vars(rep)
              }else{
                df_plotp_merged_rep <- df_plotp_smoothed

                row_var <- vars(compare_group)
              }

              # collapse
              if(collapse == TRUE){
                facet_layer <- do.call(facet_grid,modifyList(list(rows = NULL,
                                                                  cols = vars(gene_name),
                                                                  scales = "fixed"),
                                                             facet_grid_list))
              }else{
                facet_layer <- do.call(facet_grid,modifyList(list(rows = row_var,
                                                                  cols = vars(gene_name),
                                                                  scales = "fixed"),
                                                             facet_grid_list))
              }

              # ================================================================
              # plot
              # ================================================================
              if(!is.null(sample_color)){
                color_layer <- scale_color_manual(values = sample_color)
                fill_layer <- scale_fill_manual(values = sample_color)
              }else{
                color_layer <- NULL
                fill_layer <- NULL
              }

              if(mode == "codon"){
                xlabel <- "Ribosome along transcript position\n (codons / amino acids)"
              }else{
                xlabel <- "Ribosome along transcript position\n (nucleotides / nts)"
              }

              pmain <-
                ggplot(df_plotp_merged_rep) +
                # geom_ribbon(aes(x = pos,y = ratio,
                #                 ymin = ratio - mean_sd,
                #                 ymax = ratio + mean_sd,fill = compare_group),
                #             alpha = 0.5) +
                geom_plot_layer +
                theme_bw() +
                # facet_grid(rows = vars(rep),cols = vars(gene_name),scales = "fixed") +
                facet_layer +
                theme(axis.text = element_text(colour = "black"),
                      panel.grid = element_blank(),
                      strip.placement = "outside",
                      strip.text = element_text(face = "bold",size = rel(1)),
                      strip.text.x = element_text(face = "bold.italic",size = rel(1)),
                      strip.background = element_blank(),
                      strip.text.y.right = strip_text_y) +
                color_layer + fill_layer +
                xlab(xlabel) +
                ylab(ylabel) +
                xlim(c(0,cds_len))


              # ==========================================================================
              # add structure
              # ==========================================================================
              cds_rg_struc <- data.frame(x = 0,xend = cds_len,
                                         y = 1,yend = 1)

              # whether add new structure
              if(!is.null(cds_region_anno)){
                tmp_struc <- subset(cds_region_anno,gene_name == select_gene[x])
                struc_anno <- data.frame(x = as.numeric(tmp_struc$start),
                                         xend = as.numeric(tmp_struc$end),
                                         y = 1,yend = 1)

                struc_anno_layer <- ggside::geom_xsidesegment(data = struc_anno,
                                                              mapping = aes(x = x,xend = xend,y = y,yend = yend),
                                                              linewidth = cds_region_width,color = cds_anno_col)
              }else{
                struc_anno_layer <- NULL
              }

              # plot
              pmain_anno <- pmain +
                ggside::geom_xsidesegment(data = cds_rg_struc,
                                          mapping = aes(x = x,xend = xend,y = y,yend = yend),
                                          linewidth = cds_region_width,color = cds_col) +
                struc_anno_layer +
                theme(ggside.panel.background = element_blank(),
                      ggside.panel.border = element_blank(),
                      ggside.panel.scale.x = cds_region_scale_x) +
                ggside::scale_xsidey_continuous(breaks = NULL) +
                ggside::ggside(x.pos = "top",collapse = "x")

              # ==========================================================================
              # add frame annotation
              # ==========================================================================
              if(add_frame_anno == TRUE & mode == "nt"){
                frame_plot <- seq_frame_plot(cds_seq = object@CDS.sequence,
                                             trans_id = gene_df$trans_id,
                                             seq_mark = seq_mark)

                pmain_anno_tmp <- pmain_anno + xlab("") +
                  theme(plot.margin = margin(b = 0,unit = "cm"),
                        axis.title.x = element_blank())

                pmain_anno_tmp %>%
                  aplot::insert_bottom(plot = frame_plot$plot +
                                         xlab("Ribosome along transcript position\n (codons / amino acids)"),
                                       height = frame_relheight)

              }else{
                return(pmain_anno)
              }
            }) -> plist

            # combine plot list
            Reduce("+",plist)
          }
)
