globalVariables(c("RPM"))


#' Load RNA Coverage Data
#'
#' This function loads RNA coverage data from a specified file or files, processes it,
#' and returns a data frame with additional metadata including sample name, replicate, and group information.
#'
#' @param trans_density_file A character vector of file names that contain RNA transcript density data.
#' Each file should be in a tab-delimited format.
#' @param sample_name A character vector of sample names corresponding to each file in `trans_density_file`.
#' @param replicate A character vector indicating replicate information for each file in `trans_density_file`.
#' Defaults to `NULL`, in which case replicate information is set to `NA`.
#' @param group A character vector indicating the group to which each sample belongs.
#' Defaults to `NULL`, in which case group information is set to `NA`.
#'
#' @return A data frame containing the processed RNA coverage data, with columns for transcript ID (`trans_id`),
#' transcript position (`trans_pos`), reads per million (`RPM`), sample name (`sample`), replicate information (`rep`),
#' group information (`group`), and a type column (`type`) indicating the data type as RNA.
#'
#' @details The function reads in tab-delimited RNA transcript density data from the specified files, filters
#' out rows where RPM is zero, and adds metadata columns such as sample name, replicate, and group. The resulting
#' data frame is created by binding the processed data from all input files together.
#'
#'
#' @export
load_rna_coverage <- function(trans_density_file = NULL,
                              sample_name = NULL,
                              replicate = NULL,
                              group = NULL){
  # x = 1
  purrr::map_df(seq_along(trans_density_file),function(x){
    rna_tmp <- vroom::vroom(paste('2.density-data/',trans_density_file[x],sep = ''),
                            delim = "\t",show_col_types = FALSE,
                            col_names = c('gene_name',"trans_id",'trans_pos','RPM')) %>%
      dplyr::filter(RPM > 0) %>%
      dplyr::select(-gene_name) %>%
      dplyr::mutate(sample = sample_name[x],.before = trans_id) %>%
      dplyr::mutate(type = "rna",.after = RPM)

    if(!is.null(replicate)){
      rna_tmp <- rna_tmp %>%
        dplyr::mutate(rep = replicate[x],.after = sample)
    }else{
      rna_tmp <- rna_tmp %>%
        dplyr::mutate(rep = NA,.after = sample)
    }

    if(!is.null(group)){
      rna_tmp <- rna_tmp %>%
        dplyr::mutate(group = group[x],.after = rep)
    }else{
      rna_tmp <- rna_tmp %>%
        dplyr::mutate(group = NA,.after = rep)
    }

    return(rna_tmp)
  }) -> rna_df
}






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
#' @param aes_col Aesthetic var for fill color, "sample" or "type".
#' @param rna_track_df A data frame contains RNA coverage data which from load_rna_coverage function.
#' @param scale_RNA_factor Relative scale density size factor for RNA for better visualiztion. Default is 1.
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
                    aes_col = c("sample","type"),
                    rna_track_df = NULL,
                    scale_RNA_factor = 1,
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
                   aes_col = c("sample","type"),
                   rna_track_df = NULL,
                   scale_RNA_factor = 1,
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
            aes_col <- match.arg(aes_col,c("sample","type"))
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
              dplyr::mutate(type = "ribo",.after = RPM)

            # check RNA track
            if(!is.null(rna_track_df)){
              ft_rna_df <- subset(rna_track_df,trans_id == ref$tid & sample %in% unique(g_df$sample)) %>%
                dplyr::mutate(RPM = scale_RNA_factor*RPM)

              g_df <- rbind(g_df,ft_rna_df)
            }

            g_df <- g_df %>%
            dplyr::left_join(y = ref,by = c("trans_id" = "tid")) %>%
              dplyr::mutate(pos = trans_pos - utr5)

            # check mode
            if(mode == "codon"){
              g_df <- g_df %>%
                dplyr::filter(pos > 0) %>%
                dplyr::mutate(pos = ifelse(pos %% 3 == 0,pos/3,pos%/%3 + 1)) %>%
                dplyr::group_by(sample,rep,group,type,gene_name,trans_id,pos) %>%
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
                dplyr::group_by(rep,group,type,gene_name,trans_id,pos) %>%
                dplyr::summarise(RPM = mean(RPM),.groups = "drop") |>
                dplyr::rename(sample = rep)
            }

            # ==================================================================
            # main plot
            # ==================================================================
            if(!is.null(rna_track_df)){
              col_rna_layer <- do.call(geom_col,modifyList(list(data = subset(g_df,type == "rna"),
                                                                mapping = aes(x = pos,y = RPM),
                                                                width = 1),
                                                           geom_col_params))

              col_ribo_layer <- do.call(geom_col,modifyList(list(data = subset(g_df,type == "ribo"),
                                                                 mapping = aes(x = pos,y = RPM),
                                                                 width = 1),
                                                            geom_col_params))
            }else{
              col_ribo_layer <- do.call(geom_col,modifyList(list(data = g_df,
                                                                 mapping = aes(x = pos,y = RPM),
                                                                 width = 1),
                                                            geom_col_params))

              col_rna_layer <- NULL
            }

            pmian <-
              ggplot(mapping = aes(fill = get(aes_col))) +
              # geom_col(aes(x = pos,y = RPM),width = 1) +
              col_rna_layer +
              col_ribo_layer +
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
