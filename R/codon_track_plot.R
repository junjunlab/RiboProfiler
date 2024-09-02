#' Generate Codon Track Plot
#'
#' This function creates a track plot of codon-level ribosome footprint density for selected genes.
#'
#' @param codon_exp_file A vector of file paths to codon expression data files.
#' @param longest_trans_file A string specifying the path to the file containing the longest transcript information.
#' @param sample_name A vector of sample names corresponding to each codon expression file.
#' @param group_name An optional vector of group names for each sample.
#' @param select_gene A vector of gene names to be plotted.
#' @param geom_col_list A list of additional parameters to be passed to geom_col().
#'
#' @return A ggplot object representing the combined track plot for all selected genes.
#'
#' @details This function performs the following steps:
#'   1. Loads codon expression data from multiple files.
#'   2. Processes gene annotation data.
#'   3. Filters and combines data for selected genes.
#'   4. Creates individual track plots for each gene.
#'   5. Combines all plots into a single ggplot object.
#'
#' The resulting plot shows ribosome footprint density along transcript positions for each selected gene,
#' with separate tracks for each sample or group.
#'
#' @examples
#' \dontrun{codon_track_plot(
#'   codon_exp_file = c("sample1_codon_exp.txt", "sample2_codon_exp.txt"),
#'   longest_trans_file = "longest_transcripts.txt",
#'   sample_name = c("Sample1", "Sample2"),
#'   select_gene = c("GENE1", "GENE2"),
#'   geom_col_list = list(width = 0.8, position = "dodge")
#' )}
#'
#'
#' @export
codon_track_plot <- function(codon_exp_file = NULL,
                             longest_trans_file = NULL,
                             sample_name = NULL,
                             group_name = NULL,
                             select_gene = NULL,
                             geom_col_list = list()){
  # ============================================================================
  # load codon exp data
  # ============================================================================
  purrr::map_df(seq_along(codon_exp_file),function(x){
    codon_exp <- vroom::vroom(file = codon_exp_file[x],col_names = F,show_col_types = F)
    colnames(codon_exp) <- c("trans_id","pos","density")

    codon_exp$sample <- sample_name[x]

    if(!is.null(group_name)){
      tmp$group <- group_name[x]
    }

    return(codon_exp)
  }) -> df_codon_exp

  # ============================================================================
  # gene annotation
  # ============================================================================
  ganao <- read.delim(longest_trans_file,header = F)
  colnames(ganao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                       "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")
  ganao <- ganao[,c("gene_name","gene_id","trans_id")]
  ganao <- subset(ganao,gene_name %in% select_gene)

  # ============================================================================
  # generate track plot
  # ============================================================================
  df_plot <- df_codon_exp %>%
    dplyr::filter(trans_id %in% ganao$trans_id) %>%
    dplyr::left_join(y = ganao,by = "trans_id")

  # average replicates
  if(!is.null(group_name)){
    df_plot <- df_plot %>%
      dplyr::group_by(group,trans_id,pos) %>%
      dplyr::summarise(density = mean(density)) %>%
      dplyr::rename(sample = group)
  }

  # loop plot
  lapply(seq_along(select_gene),function(x){
    if(x != length(select_gene)){
      strip_text_y <- element_blank()
    }else{
      strip_text_y <- element_text(angle = 0)
    }


    geom_col_layer <- do.call(geom_col,
                              modifyList(list(mapping = aes(x = pos,y = density,
                                                            fill = sample,color = sample),
                                              show.legend = FALSE),
                                         geom_col_list))

    # plot
    ggplot(df_plot %>% dplyr::filter(gene_name == select_gene[x])) +
      # geom_col(aes(x = pos,y = density)) +
      geom_col_layer +
      theme_bw() +
      facet_grid(rows = vars(sample),cols = vars(gene_name),scales = "fixed") +
      theme(axis.text = element_text(colour = "black"),
            panel.grid = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(face = "bold",size = rel(1)),
            strip.background = element_blank(),
            strip.text.y.right = strip_text_y) +
      xlab("Ribosome along transcript position\n (codon)") +
      ylab("Ribosome footprint density")

  }) -> plist

  # combine plot list
  Reduce("+",plist)
}
