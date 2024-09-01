#' Calculate Relative Distance of Motifs within Transcripts
#'
#' This function calculates the relative distances of specified motifs within transcripts
#' from normalized data files, filtered by minimum expression counts, and outputs the results.
#' It handles multiple samples and uses external Python scripts for detailed calculations.
#'
#' @param amino_file Character, path to the amino acid sequence file, default is NULL.
#' @param longest_trans_file Character, path to the file with longest transcripts information, default is NULL.
#' @param normed_file Data frame, normalized data file containing sample, trans_id, and counts columns, default is NULL.
#' @param min_counts Integer, minimum count threshold for filtering transcripts, default is 64.
#' @param motif Character, the specific motif to search for, default is NULL.
#' @param upstream Integer, upstream distance to consider relative to motif position, default is -50.
#' @param downstream Integer, downstream distance to consider relative to motif position, default is 50.
#' @return Invisible NULL, as the function primarily writes output to files.
#' @examples
#' \dontrun{rel_dist_motif(amino_file = "amino_seq.txt",
#'                longest_trans_file = "longest_transcripts.txt",
#'                normed_file = data.frame(sample = c("sample1", "sample2"),
#'                                         trans_id = c("trans1", "trans2"),
#'                                         counts = c(100, 150)),
#'                min_counts = 64,
#'                motif = "PKP",
#'                upstream = -50,
#'                downstream = 50)}
#'
#' @export
rel_dist_motif <- function(amino_file = NULL,
                           longest_trans_file = NULL,
                           normed_file = NULL,
                           min_counts = 64,
                           motif = NULL,
                           upstream = -50,
                           downstream = 50){
  # =====================================================================================
  # output normalized files
  # =====================================================================================
  dir.create("relative_distance_motif",showWarnings = FALSE)

  sp <- unique(normed_file$sample)

  # loop output
  # x = 1
  lapply(seq_along(sp),function(x){
    tmp <- subset(normed_file,sample == sp[x])

    # filter low counts
    total_exp <- tmp %>%
      dplyr::group_by(trans_id) %>%
      dplyr::summarise(total_counts = sum(counts)) %>%
      dplyr::filter(total_counts >= min_counts)

    filtered_df <- tmp %>%
      dplyr::filter(trans_id %in% total_exp$trans_id)

    # output
    vroom::vroom_write(x = filtered_df,
                       file = paste("relative_distance_motif/",sp[x],"_normed.txt",sep = ""),
                       col_names = F)

    return(NULL)
  }) -> none


  # =====================================================================================
  # calculation
  # =====================================================================================

  # run code
  pyscript.path = system.file("extdata", "relDistTripeptideMotif.py", package = "RiboProfiler")
  reticulate::source_python(pyscript.path)

  # loop calculation
  # x = 1
  lapply(seq_along(sp),function(x){
    normed_file <- paste("relative_distance_motif/",sp[x],"_normed.txt",sep = "")
    output_file <- paste("relative_distance_motif/",sp[x],"_",motif,"_motif_occupancy.txt",sep = "")

    # run
    suppressMessages(
      reticulate::py$relDistTripeptideMotif(amino_file = amino_file,
                                            longest_trans_file = longest_trans_file,
                                            normed_file = normed_file,
                                            motif = motif,
                                            output_file = output_file,
                                            upstream = as.integer(upstream),
                                            downstream = as.integer(downstream))
    )

    message(paste("#",normed_file,"has been processed!"))
  }) -> none

}






#' Plot Relative Motif Occupancy from Motif Density Files
#'
#' This function takes the output from motif density calculations and generates
#' a plot of relative footprint density across positions relative to a specified motif.
#' It supports handling multiple samples and grouping them if group information is provided.
#'
#' @param motif_occupancy_file Character vector, paths to the motif occupancy files.
#' @param sample_name Character vector, names of the samples corresponding to the files.
#' @param group_name Optional character vector, names of groups for each sample.
#'        If provided, averages density across samples in the same group.
#' @param motif Character, the specific motif being analyzed.
#' @param site Character, label for the site being considered in the motif analysis,
#'        default is "P" (e.g., for the Peptidyl site in ribosome profiling).
#' @return An ggplot object showing the relative motif occupancy plot.
#' @examples
#' \dontrun{rel_dist_motif_plot(motif_occupancy_file = c("sample1.txt", "sample2.txt"),
#'                     sample_name = c("sample1", "sample2"),
#'                     group_name = c("group1", "group1"),
#'                     motif = "PKP",
#'                     site = "P")}
#'
#' @export
rel_dist_motif_plot <- function(motif_occupancy_file = NULL,
                                sample_name = NULL,
                                group_name = NULL,
                                motif = NULL,
                                site = "P"){
  # ============================================================================
  # load files
  purrr::map_df(seq_along(motif_occupancy_file),function(x){
    rel_motif <- read.delim(motif_occupancy_file[x],header = F)
    colnames(rel_motif) <- c("pos","density")
    rel_motif$sample <- sample_name[x]

    if(!is.null(group_name)){
      tmp$group <- group_name[x]
    }

    return(rel_motif)
  }) -> motif_df

  # ============================================================================
  # deal with replicates
  if(!is.null(group_name)){
    df_plot <- motif_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(density = mean(density)) %>%
      dplyr::rename(sample = group)
  }else{
    df_plot <- motif_df
  }

  # ============================================================================
  # plot
  ggplot(motif_df,aes(x = pos,y = density,color = sample)) +
    geom_line(linewidth = 0.5) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.text = element_text(face = "bold",size = rel(1))) +
    xlab(paste("Distance of", site, "site to",motif,"(nt)")) +
    ylab("Relative footprint density (AU)") +
    scale_color_brewer(palette = "Set1")
}
