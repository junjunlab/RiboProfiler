globalVariables(c("Offsets", "abbreviation", "amino", "bamFiles", "bamLegends", "cdsft", "codon", "norm_exp",
                  "readLengths", "total_counts", "trans_pos", "anno", "relsp","..."))

#' Calculate Codon Occupancy
#'
#' This function calculates codon occupancy based on ribosome profiling data and transcript information.
#'
#' @param object ribosomeObj object.
#' @param min_counts An integer specifying the minimum number of counts required for a transcript to be included
#'        in the analysis. Default is 32.
#' @param norm_type The nomalization methods for ribosome density. "average" is calculated by
#' the count at each position divided by mean density across cds region. "rpm"
#' is calculated by the count at each position divided by the total counts and multiplied with 10^6.
#' Default is "average".
#' @param upstream_codon_exclude An integer specifying the number of codons to exclude from the start of the coding sequence.
#'        Default is 5.
#' @param downstream_codon_exclude An integer specifying the number of codons to exclude from the end of the coding sequence.
#'        Default is 0.
#' @param ... Useless args.
#'
#' @details This function performs the following steps:
#'   1. Filters reads based on minimum count threshold.
#'   2. Calculates codon density for each transcript.
#'   3. Outputs codon position expression data.
#'   4. Calls a Python script to calculate codon occupancy.
#'
#' @return This function does not return a value, but creates output files in a 'codon_occupancy' directory:
#'   - A text file for each sample containing codon position expression data.
#'   - A text file for each sample containing codon occupancy data.
#'
#'
#' @importFrom stats na.omit
#'
#' @export
setGeneric("codon_occupancy",
           function(object,
                    min_counts = 32,
                    norm_type = c("average","rpm"),
                    upstream_codon_exclude = 0,
                    downstream_codon_exclude = 0, ...) standardGeneric("codon_occupancy"))





#' method for codon_occupancy
#'
#' @rdname codon_occupancy
#' @exportMethod codon_occupancy
setMethod("codon_occupancy",
          signature(object = "ribosomeObj"),
          function(object,
                   min_counts = 32,
                   norm_type = c("average","rpm"),
                   upstream_codon_exclude = 0,
                   downstream_codon_exclude = 0,...){
            norm_type <- match.arg(norm_type,c("average","rpm"))
            # ============================================================================
            # gene annotation
            # ============================================================================
            dir.create("codon_occupancy",showWarnings = FALSE)

            ganao <- object@longest.annotation
            colnames(ganao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                                 "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

            # loop for each sample
            normed_file <- object@normalized.data
            sp <- unique(normed_file$sample)

            # x = 1
            lapply(seq_along(sp),function(x){
              # ============================================================================
              # filter reads by counts
              # ============================================================================
              inputfile <- normed_file %>%
                dplyr::filter(sample == sp[x])

              total_exp <- inputfile %>%
                dplyr::left_join(y = ganao,by = "trans_id") %>%
                dplyr::filter(trans_pos >= `5UTR_length` & trans_pos <= `5UTR_length` + CDS_length) %>%
                dplyr::group_by(trans_id) %>%
                dplyr::summarise(total_counts = sum(counts)) %>%
                dplyr::filter(total_counts >= min_counts)

              ganao <- ganao[,c("trans_id","5UTR_length","CDS_length","3UTR_length")] %>%
                dplyr::filter(trans_id %in% total_exp$trans_id)

              # ============================================================================
              # get codon density
              # ============================================================================

              qc_df_anno <- inputfile %>%
                dplyr::filter(trans_id %in% total_exp$trans_id) %>%
                dplyr::left_join(y = ganao,by = "trans_id") %>%
                dplyr::mutate(cdsft = dplyr::if_else(CDS_length%%3 == 0,1,0)) %>%
                # filter cds length %%3 == 0
                dplyr::filter(cdsft == 1) %>%
                dplyr::mutate(trans_pos = trans_pos - `5UTR_length`) %>%
                # fiter trans_pos in cds region
                dplyr::filter(trans_pos > upstream_codon_exclude*3  & trans_pos <= CDS_length - downstream_codon_exclude*3) %>%
                # transpos trans_pos into codon pos
                dplyr::mutate(codon_pos = dplyr::if_else(trans_pos%%3 == 0, trans_pos/3,trans_pos%/%3 + 1)) %>%
                dplyr::group_by(trans_id,codon_pos) %>%
                dplyr::summarise(codon_exp = mean(!!rlang::sym(norm_type))) %>%
                na.omit()

              # output
              fname = paste("codon_occupancy/",sp[x],"_codon_pos_exp.txt",sep = "")
              vroom::vroom_write(x = qc_df_anno,file = fname,col_names = F,quote = "none")

              # ============================================================================
              # get codon seq and average density
              # ============================================================================
              # run code
              pyscript.path = system.file("extdata", "CodonOccupancy.py", package = "RiboProfiler")
              reticulate::source_python(pyscript.path)

              if(length(object@CDS.sequence) != 0){
                suppressMessages(
                  reticulate::py$codonOccupancy(cds_fasta_file = object@CDS.sequence,
                                                codon_pos_exp_file = fname,
                                                output_file = paste("codon_occupancy/",sp[x],"_codon_occupancy.txt",sep = ""))
                )
              }else{
                message("Please run fetch_sequence first to get cds fasta file!")
              }

            })
          }
)








#' Plot Codon Occupancy
#'
#' This function creates a plot of codon occupancy based on the provided data.
#'
#' @param codon_occupancy_file A vector of file paths to codon occupancy data files.
#' @param sample_name A vector of sample names corresponding to each file.
#' @param group_name An optional vector of group names for each sample.
#' @param compare_var A vector of two variable names to compare.
#' @param facet A logical value indicating whether to use faceting in the plot. Default is FALSE.
#' @param codon_type A character string specifying whether to plot by codon or amino acid.
#'        Options are "codon" or "amino". Default is c("codon", "amino").
#' @param rm_stopcodon Whether remove stop codons. Default is FALSE.
#' @param col The color fill types for plot. Default is c("normal","chemistry"),
#'
#' @return A ggplot object representing the codon occupancy plot.
#'
#' @details This function reads codon occupancy data from files, processes it, and creates
#'          a bar plot showing the difference in codon occupancy between two conditions.
#'          It can plot data by codon or by amino acid, and offers options for faceting.
#'
#' @examples
#' \dontrun{# Basic usage
#' codon_occupancy_plot(
#'   codon_occupancy_file = c("sample1.txt", "sample2.txt"),
#'   sample_name = c("Sample1", "Sample2"),
#'   compare_var = c("Sample1", "Sample2")
#' )
#'
#' # With grouping and faceting
#' codon_occupancy_plot(
#'   codon_occupancy_file = c("sample1.txt", "sample2.txt", "sample3.txt", "sample4.txt"),
#'   sample_name = c("Sample1", "Sample2", "Sample3", "Sample4"),
#'   group_name = c("GroupA", "GroupA", "GroupB", "GroupB"),
#'   compare_var = c("GroupA", "GroupB"),
#'   facet = TRUE,
#'   codon_type = "codon"
#' )}
#'
#' @importFrom ggh4x weave_factors
#'
#'
#' @export
codon_occupancy_plot <- function(codon_occupancy_file = NULL,
                                 sample_name = NULL,
                                 group_name = NULL,
                                 compare_var = NULL,
                                 rm_stopcodon = FALSE,
                                 col = c("normal","chemistry"),
                                 facet = FALSE,
                                 codon_type = c("codon","amino")){
  col <- match.arg(col,c("normal","chemistry"))
  codon_type <- match.arg(codon_type,c("codon","amino"))
  # ============================================================================
  # loop load file
  # ============================================================================
  purrr::map_df(seq_along(codon_occupancy_file),function(x){
    tmp <- read.delim(codon_occupancy_file[x],header = F)
    colnames(tmp) <- c("codon","amino","abbreviation","density")

    # whether remove stop codon
    if(rm_stopcodon == TRUE){
      tmp <- subset(tmp,amino != "Stop")
    }

    tmp$sample <- sample_name[x]

    if(!is.null(group_name)){
      tmp$group_name <- group_name[x]
    }

    return(tmp)
  }) -> codon_oc

  # ============================================================================
  # check type
  # ============================================================================
  ylab <- paste("Difference of codon occupancy","\n",
                "(",compare_var[1]," VS ",compare_var[2],")",sep = "")

  if(codon_type == "amino"){
    # prepare amino acids annotations
    chemistry = data.frame(
      abbreviation = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
      anno = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
      stringsAsFactors = F)

    if(!is.null(group_name)){
      df_codon <- codon_oc %>%
        group_by(group_name,abbreviation) %>%
        summarise(density = mean(density)) %>%
        tidyr::spread(group_name,density)

    }else{
      df_codon <- codon_oc %>%
        group_by(sample,abbreviation) %>%
        summarise(density = mean(density)) %>%
        tidyr::spread(sample,density)
    }

    df_codon <- df_codon %>%
      dplyr::left_join(y = chemistry,by = "abbreviation")

    df_codon <- df_codon %>%
      dplyr::mutate(ratio = log2(.data[[compare_var[1]]]/.data[[compare_var[2]]])) %>%
      dplyr::arrange(dplyr::desc(ratio))

    if(col == "chemistry"){
      df_codon <- df_codon %>% na.omit()
      col_layer <- geom_col(aes(x = abbreviation,y = ratio,fill = anno),
                            show.legend = T,width = 0.5)

      amino_cols <- c('#109648','#5E239D','#255C99','#D62839','#221E22')
      names(amino_cols) <- c("Polar","Neutral","Basic","Acidic","Hydrophobic")
      color_layer <- scale_fill_manual(values = amino_cols,name = "Amino acids types")
    }else{
      col_layer <- geom_col(aes(x = abbreviation,y = ratio,fill = abbreviation),
                            show.legend = F,width = 0.5)

      color_layer <- scale_color_brewer()
    }

    # order
    df_codon$abbreviation <- factor(df_codon$abbreviation,levels = df_codon$abbreviation)

    # plot
    ggplot(df_codon) +
      # geom_col(aes(x = abbreviation,y = ratio,fill = abbreviation),
      #          show.legend = F,width = 0.5) +
      col_layer +
      geom_hline(yintercept = 0,lty = "solid",color = "black",linewidth = 0.75) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text = element_text(color = "black")) +
      ylab(ylab) +
      xlab("Codons / Amino acids") +
      color_layer
  }else{
    if(!is.null(group_name)){
      df_codon <- codon_oc %>%
        group_by(group_name,amino,codon) %>%
        summarise(density = mean(density)) %>%
        tidyr::spread(group_name,density)
    }else{
      df_codon <- codon_oc %>%
        tidyr::spread(sample,density)
    }

    df_codon <- df_codon %>%
      dplyr::mutate(ratio = log2(.data[[compare_var[1]]]/.data[[compare_var[2]]])) %>%
      dplyr::arrange(dplyr::desc(ratio))

    # order
    df_codon$codon <- factor(df_codon$codon,levels = df_codon$codon)

    # plot
    if(facet == FALSE){
      ggplot(df_codon) +
        geom_col(aes(x = ggh4x::weave_factors(codon,amino),y = ratio,fill = codon),
                 width = 0.5,show.legend = F) +
        geom_hline(yintercept = 0,lty = "solid",color = "black",linewidth = 0.75) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(color = "black"),
              ggh4x.axis.nesttext.x = element_text(angle = 90,hjust = 1),
              axis.text.x = element_text(angle = 90,vjust = 0.5)) +
        guides(x = "axis_nested") +
        ylab(ylab) +
        xlab("Codons / Amino acids")
    }else{
      df_codon$amino <- factor(df_codon$amino,levels = unique(df_codon$amino))

      ggplot(df_codon) +
        geom_col(aes(x = codon,y = ratio,fill = amino),
                 width = 0.5,show.legend = F) +
        geom_hline(yintercept = 0,lty = "solid",color = "black",linewidth = 0.75) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              strip.text = element_text(face = "bold",size = rel(1)),
              axis.text = element_text(color = "black"),
              ggh4x.axis.nesttext.x = element_text(angle = 90),
              axis.text.x = element_text(angle = 90,vjust = 0.5)) +
        facet_grid(~amino,scales = "free_x",space = "free_x") +
        ylab(ylab) +
        xlab("Codons / Amino acids")
    }

  }
}

