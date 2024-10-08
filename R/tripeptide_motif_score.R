globalVariables(c("motif", "score"))

#' Calculate Peptide Motif Scores
#'
#' This function calculates peptide motif scores based on amino acid sequences and codon expression data.
#' It uses a Python script to perform the calculations and saves the results for each sample.
#'
#' @param object ribosomeObj object.
#' @param min_counts Minimum counts to filter gene, default 64.
#' @param norm_type The nomalization methods for ribosome density. "average" is calculated by
#' the count at each position divided by mean density across cds region. "rpm"
#' is calculated by the count at each position divided by the total counts and multiplied with 10^6.
#' Default is "average".
#' @param average_normalization Whether do average normalization, yes or no.
#' @param occurrence_threshold An integer specifying the minimum occurrence threshold for peptide motifs.
#' @param ... Useless args.
#'
#' @return This function does not return a value, but creates output files in a 'peptide_motif' directory.
#'         Each output file is named '[sample_name]_tripeptide_occupancy.txt'.
#'
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
setGeneric("peptide_motif_score",
           function(object,
                    min_counts = 64,
                    norm_type = c("average","rpm"),
                    average_normalization = c("yes","no"),
                    occurrence_threshold = 100, ...) standardGeneric("peptide_motif_score"))




#' method for peptide_motif_score
#'
#' @rdname peptide_motif_score
#' @exportMethod peptide_motif_score
setMethod("peptide_motif_score",
          signature(object = "ribosomeObj"),
          function(object,
                   min_counts = 64,
                   norm_type = c("average","rpm"),
                   average_normalization = c("yes","no"),
                   occurrence_threshold = 100,...){
            average_normalization <- match.arg(average_normalization,c("yes","no"))
            norm_type <- match.arg(norm_type,c("average","rpm"))
            norm_type <- ifelse(norm_type == "average","rpm","average")
            # =====================================================================================
            # calculation
            # =====================================================================================
            dir.create("peptide_motif",showWarnings = FALSE)

            # source code
            pyscript.path = system.file("extdata", "peptideMotifScore.py", package = "RiboProfiler")
            reticulate::source_python(pyscript.path)

            normed_exp_file <- object@normalized.data
            sp <- unique(normed_exp_file$sample)
            # x = 1
            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(normed_exp_file,sample == sp[x]) %>%
                dplyr::select(-dplyr::all_of(norm_type))

              vroom::vroom_write(tmp,col_names = F,
                                 file = paste("peptide_motif/",sp[x],"_normed.txt",sep = ""))

              # calculation
              # check amino acid file
              if(length(object@amino.acid.sequence) != 0){
                suppressMessages(
                  reticulate::py$peptideMotifScore(amino_file = object@amino.acid.sequence,
                                                   longest_trans_file = object@longest.anno.file,
                                                   normed_file = paste("peptide_motif/",sp[x],"_normed.txt",sep = ""),
                                                   min_counts = as.integer(min_counts),
                                                   average_normalization = average_normalization,
                                                   output_file = paste("peptide_motif/",sp[x],
                                                                       "_tripeptide_occupancy.txt",sep = ""),
                                                   occurrence_threshold = as.integer(occurrence_threshold))
                )

                message(paste(sp[x],"has been processed!"))
              }else{
                message("Please run fetch_sequence first to get amino acids file!")
              }

              return(NULL)
            }) -> tmp_normed

          }
)






#' Peptide Motif Scoring Function
#'
#' This function calculates the motif scores for peptides based on provided files.
#'
#' @param object ribosomeObj object.
#' @param norm_type A character string specifying the normalization type. It can be one of
#'   "rpm" or "average".
#'   Default is "normed_count".
#' @param window An integer specifying the window size for motif scanning.
#'   Default is 50.
#' @param occurrence_threshold An integer specifying the threshold for motif occurrences.
#'   Default is 100.
#' @param ... Useless args.
#'
#' @return The function creates output files with motif scores and prints a message for each processed sample.
#'
#' @examples
#' \dontrun{# Example usage:
#' peptide_motif_score2(object = ribosomeObj)}
#'
#' @export
setGeneric("peptide_motif_score2",
           function(object,
                    norm_type = c("rpm","average"),
                    window = 50,
                    occurrence_threshold = 100, ...) standardGeneric("peptide_motif_score2"))





#' method for peptide_motif_score2
#'
#' @rdname peptide_motif_score2
#' @exportMethod peptide_motif_score2
setMethod("peptide_motif_score2",
          signature(object = "ribosomeObj"),
          function(object,
                   norm_type = c("rpm","average"),
                   window = 50,
                   occurrence_threshold = 100,...){
            norm_type <- match.arg(norm_type,c("rpm","average"))
            norm_type <- ifelse(norm_type == "average","rpm","average")
            # =====================================================================================
            # output
            # =====================================================================================
            dir.create("averaged_peptide_motif_data",showWarnings = F)

            # source code
            pyscript.path = system.file("extdata", "averageMotifPauseScore.py", package = "RiboProfiler")
            reticulate::source_python(pyscript.path)

            normed_exp_file <- object@normalized.data
            sp <- unique(normed_exp_file$sample)
            # x = 1
            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(normed_exp_file,sample == sp[x]) %>%
                dplyr::select(-dplyr::all_of(norm_type))


              vroom::vroom_write(tmp,col_names = F,
                                 file = paste("averaged_peptide_motif_data/",sp[x],"_normed.txt",sep = ""))

              # calculation
              # check amino acid file
              if(length(object@amino.acid.sequence) != 0){
                suppressMessages(
                  reticulate::py$averageMotifPauseScore(amino_file = object@amino.acid.sequence,
                                                        longest_trans_file = object@longest.anno.file,
                                                        input_file = paste("averaged_peptide_motif_data/",sp[x],"_normed.txt",sep = ""),
                                                        output_file = paste("averaged_peptide_motif_data/",sp[x],
                                                                            "_tripeptide_occupancy.txt",sep = ""),
                                                        norm_type = norm_type,
                                                        window = as.integer(window),
                                                        occurrence_threshold = as.integer(occurrence_threshold))
                )

                message(paste(sp[x],"has been processed!"))
              }else{
                message("Please run fetch_sequence first to get amino acids file!")
              }

              return(NULL)
            }) -> tmp_normed
          }
)



#' Create a Scatter Plot of Tripeptide Occupancy Scores
#'
#' This function generates a scatter plot comparing tripeptide occupancy scores between two samples or groups.
#'
#' @param occupancy_file A vector of file paths to tripeptide occupancy score files.
#' @param sample_name A vector of sample names corresponding to each occupancy file.
#' @param group_name An optional vector of group names for each sample, it is useful
#' to deal with replicates by using mean method.
#' @param mark_motif A vector of motifs to highlight in the plot.
#' @param mark_color The color to use for highlighted motifs. Default is "orange".
#' @param pcol The point colors. Default is "grey".
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
                                 sample_name = NULL,
                                 group_name = NULL,
                                 mark_motif = NULL,
                                 pcol = "grey",
                                 mark_color = "orange"){
  purrr::map_df(seq_along(occupancy_file),function(x){
    tmp <- read.delim(occupancy_file[x],header = F)
    colnames(tmp) <- c("motif","score")
    tmp$score <- as.numeric(tmp$score)

    tmp$sample <- sample_name[x]

    if(!is.null(group_name)){
      tmp$group <- group_name[x]
    }

    return(tmp)
  }) -> df_plot

  # get mean score for different groups
  if(!is.null(group_name)){
    var <- unique(group_name)
    df_plot_wide <- df_plot %>%
      dplyr::group_by(group,motif) %>%
      dplyr::summarise(score = mean(score)) %>%
      tidyr::spread(group,score) %>%
      na.omit()
  }else{
    var <- sample_name
    df_plot_wide <- df_plot %>%
      tidyr::spread(sample,score) %>%
      na.omit()
  }

  mark <- subset(df_plot_wide,motif %in% mark_motif)

  # plot
  ggplot(df_plot_wide) +
    geom_point(aes(x = .data[[var[1]]],y = .data[[var[2]]]),
               color = pcol) +
    geom_abline(slope = 1,lty = "dashed",linewidth = 0.75) +
    ggrepel::geom_text_repel(data = mark,
                             aes(x = .data[[var[1]]],y = .data[[var[2]]],
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






#' Plot Tri-amino Acid Motifs
#'
#' This function takes motif occupancy data and generates a plot using the ggseqlogo
#' package to visualize the frequency or probability of amino acid motifs. It handles
#' multiple samples and can average scores across specified groups.
#'
#' @param occupancy_file Vector of strings; paths to the occupancy files containing
#'        motifs and their scores.
#' @param sample_name Vector of strings; names of the samples corresponding to each
#'        file. Used for labeling in the absence of a group name.
#' @param group_name Optional; vector of strings providing group names for each sample.
#'        If provided, scores are averaged across groups.
#' @param top_motif Integer; number of top motifs to display based on their ratio.
#' @return A ggplot object displaying the sequence logo of top motifs based on the
#'         calculated ratio.
#' @examples
#' \dontrun{triAmino_motif_plot(occupancy_file = c("path/to/file1.txt", "path/to/file2.txt"),
#'                     sample_name = c("sample1", "sample2"),
#'                     group_name = c("group1", "group1"),
#'                     top_motif = 5)}
#'
#' @importFrom ggseqlogo ggseqlogo
#'
#' @export
triAmino_motif_plot <- function(occupancy_file = NULL,
                                sample_name = NULL,
                                group_name = NULL,
                                top_motif = NULL){
  # x = 1
  purrr::map_df(seq_along(occupancy_file),function(x){
    tmp <- read.delim(occupancy_file[x],header = F)
    colnames(tmp) <- c("motif","score")

    tmp$sample <- sample_name[x]

    if(!is.null(group_name)){
      tmp$group <- group_name[x]
    }

    return(tmp)
  }) -> df_plot

  # get mean score for different groups
  if(!is.null(group_name)){
    var <- group_name
    df_plot_wide <- df_plot %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(score = mean(score)) %>%
      tidyr::spread(group,score) %>%
      na.omit()
  }else{
    var <- sample_name
    df_plot_wide <- df_plot %>%
      tidyr::spread(sample,score) %>%
      na.omit()
  }

  df_plot_wide <- df_plot_wide %>%
    dplyr::mutate(ratio = .data[[var[2]]]/.data[[var[1]]]) %>%
    dplyr::arrange(dplyr::desc(ratio)) %>%
    head(n = top_motif)

  # plot
  ggseqlogo::ggseqlogo(data = df_plot_wide$motif,method="prob",col_scheme = "chemistry") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black")) +
    scale_x_continuous(labels = c("E","P","A"),breaks = c(1,2,3))

}

