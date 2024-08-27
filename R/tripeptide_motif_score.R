globalVariables(c("motif", "score"))

#' Calculate Peptide Motif Scores
#'
#' This function calculates peptide motif scores based on amino acid sequences and codon expression data.
#' It uses a Python script to perform the calculations and saves the results for each sample.
#'
#' @param amino_file A string specifying the path to the file containing amino acid sequences.
#' @param codon_exp_file A vector of file paths to codon expression data files.
#' @param sample_name A vector of sample names corresponding to each codon expression file.
#' @param occurrence_threshold An integer specifying the minimum occurrence threshold for peptide motifs.
#'        Default is 50.
#'
#' @return This function does not return a value, but creates output files in a 'peptide_motif' directory.
#'         Each output file is named '[sample_name]_tripeptide_occupancy.txt'.
#'
#' @examples
#' \dontrun{peptide_motif_score(
#'   amino_file = "path/to/amino_sequences.txt",
#'   codon_exp_file = c("sample1_codon_exp.txt", "sample2_codon_exp.txt"),
#'   sample_name = c("Sample1", "Sample2"),
#'   occurrence_threshold = 50
#' )}
#'
#' @details This function performs the following steps:
#'   1. Creates a 'peptide_motif' directory if it doesn't exist.
#'   2. Sources a Python script 'peptideMotifScore.py' from the RiboProfiler package.
#'   3. For each codon expression file:
#'      - Calls the Python function to calculate peptide motif scores.
#'      - Saves the results in the 'peptide_motif' directory.
#'      - Prints a message indicating the file has been processed.
#'
#' @export
peptide_motif_score <- function(amino_file = NULL,
                                codon_exp_file = NULL,
                                sample_name = NULL,
                                occurrence_threshold = 50){
  # =====================================================================================
  # calculation
  # =====================================================================================
  dir.create("peptide_motif",showWarnings = FALSE)

  # run code
  pyscript.path = system.file("extdata", "peptideMotifScore.py", package = "RiboProfiler")
  reticulate::source_python(pyscript.path)

  # loop calculation
  # x = 2
  lapply(seq_along(codon_exp_file),function(x){
    suppressMessages(
      reticulate::py$peptideMotifScore(amino_file = amino_file,
                                       codon_exp_file = codon_exp_file[x],
                                       output_file = paste("peptide_motif/",sample_name[x],
                                                           "_tripeptide_occupancy.txt",sep = ""),
                                       occurrence_threshold = occurrence_threshold)
    )

    message(paste(codon_exp_file[x],"has been processed!"))
  })

}






#' Create a Scatter Plot of Tripeptide Occupancy Scores
#'
#' This function generates a scatter plot comparing tripeptide occupancy scores between two samples or groups.
#'
#' @param occupancy_file A vector of file paths to tripeptide occupancy score files.
#' @param sampe_name A vector of sample names corresponding to each occupancy file.
#' @param group_name An optional vector of group names for each sample.
#' @param mark_motif A vector of motifs to highlight in the plot.
#' @param mark_color The color to use for highlighted motifs. Default is "orange".
#'
#' @return A ggplot object representing the scatter plot of tripeptide occupancy scores.
#'
#' @details This function performs the following steps:
#'   1. Reads and combines tripeptide occupancy data from multiple files.
#'   2. Calculates mean scores if group names are provided.
#'   3. Creates a scatter plot comparing scores between two samples or groups.
#'   4. Highlights specified motifs with labels and different colors.
#'
#' @examples
#' \dontrun{triAmino_scater_plot(
#'   occupancy_file = c("sample1_occupancy.txt", "sample2_occupancy.txt"),
#'   sampe_name = c("Sample1", "Sample2"),
#'   mark_motif = c("PPP", "GGG"),
#'   mark_color = "red"
#' )}
#'
#' @importFrom ggrepel geom_text_repel
#'
#' @export
triAmino_scater_plot <- function(occupancy_file = NULL,
                                 sampe_name = NULL,
                                 group_name = NULL,
                                 mark_motif = NULL,
                                 mark_color = "orange"){
  purrr::map_df(seq_along(occupancy_file),function(x){
    tmp <- read.delim(occupancy_file[x],header = F)
    colnames(tmp) <- c("motif","score")

    tmp$sample <- sampe_name[x]

    if(!is.null(group_name)){
      tmp$group <- group_name[x]
    }

    return(tmp)
  }) -> df_plot

  # get mean score for different groups
  if(!is.null(group_name)){
    var <- group_name
    df_plot <- df_plot %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(score = mean(score)) %>%
      tidyr::spread(group,score) %>%
      na.omit()
  }else{
    var <- sampe_name
    df_plot_wide <- df_plot %>%
      tidyr::spread(sample,score) %>%
      na.omit()
  }

  mark <- subset(df_plot_wide,motif %in% mark_motif)

  # plot
  ggplot(df_plot_wide) +
    geom_point(aes(x = .data[[var[1]]],y = .data[[var[2]]]),
               color = "grey") +
    geom_abline(slope = 1,lty = "dashed",linewidth = 0.75) +
    ggrepel::geom_text_repel(data = mark,
                             aes(x = .data[["WT"]],y = .data[["eIF5Ad"]],
                                 label = motif),
                             color = mark_color,
                             max.overlaps = Inf,
                             direction = "both",
                             segment.size = 0.5,
                             force = 50, max.iter = 3e3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(colour = "black")) +
    coord_fixed(ratio = 1) +
    xlim(c(0,max(df_plot$score))) +
    ylim(c(0,max(df_plot$score))) +
    labs(title = "tripeptide pausing score")

}
