globalVariables(c("abs_pos", "cds", "count", "pos", "sum_density" ,
                  "sum_pi", "tid", "utr5", "wi", "x", "y"))

#' Calculate Polarity Scores
#'
#' This function calculates polarity scores for genes based on RNA-seq data. The
#' function reads gene annotation, GTF file, SAM file, density file, and optionally
#' a grouping factor. It then iterates over the provided SAM files, calculates counts
#' for each gene in the CDS region, filters the counts to those above a threshold,
#' and calculates the polarity score (PS) for each gene.
#'
#' @param gene_anno_file A character string specifying the path to the gene annotation
#'   file. Default is NULL.
#' @param gtf_file A character string specifying the path to the GTF file. Default is
#'   NULL.
#' @param input_file A character vector specifying the paths to the SAM/BAM files. Default is
#'   NULL.
#' @param density_file A character vector specifying the paths to the density files.
#'   Default is NULL.
#' @param minCounts A numeric value specifying the minimum counts threshold for
#'   filtering the counts. Default is 64.
#' @param upstreamNt A numeric value specifying the number of nucleotides upstream
#'   of the UTR5 region to include. Default is 15.
#' @param downstreamNt A numeric value specifying the number of nucleotides downstream
#'   of the CDS region to include. Default is 15.
#' @param group A character vector specifying a grouping factor. Default is NULL.
#' @return A data frame containing the polarity scores.
#' @importFrom Rsubread featureCounts
#' @importFrom dplyr filter left_join mutate group_by summarise
#' @importFrom purrr map_df
#' @importFrom vroom vroom
#' @import stringr
#' @export
calculatePolarity <- function(gene_anno_file = NULL,
                              gtf_file = NULL,
                              input_file = NULL,
                              density_file = NULL,
                              minCounts = 64,
                              upstreamNt = 15,
                              downstreamNt = 15,
                              group = NULL){
  # ========================================================================================
  # load gene annotation
  # ========================================================================================
  gene_anno <- read.delim("longest_info.txt",header = F)

  # add colnames
  colnames(gene_anno) <- c("id","gene_name","gdi","tid","chr","strand","cds_rg","exon_rg","utr5","cds","utr3")

  ganno <- gene_anno[,c("tid","utr5","cds","utr3")]

  # ========================================================================================
  # define internal function
  # ========================================================================================
  calc_PS_fun <- function(gtf_file = NULL,
                          input_file = NULL,
                          density_file = NULL,
                          minCounts = 64,
                          upstreamNt = 15,
                          downstreamNt = 15,
                          group = NULL){
    # ========================================================================================
    # calculate counts for each gene
    # ========================================================================================
    exp <- Rsubread::featureCounts(files = input_file,
                                   annot.ext = gtf_file,
                                   isGTFAnnotationFile = T,
                                   GTF.featureType = "CDS",
                                   GTF.attrType = "gene_id",
                                   GTF.attrType.extra = c("gene_name"),
                                   isPairedEnd = FALSE,
                                   nthreads = parallel::detectCores())

    # filter counts threshold in CDS for each gene
    select_gene <- data.frame(count = as.vector(exp$counts),exp$annotation) |>
      subset(count > minCounts)

    # ========================================================================================
    # polarity calculation
    # ========================================================================================
    # load density files
    dt <- vroom::vroom(file = density_file,
                       delim = "\t",col_names = FALSE,show_col_types = FALSE)

    # add colnames
    colnames(dt) <- c("gene_name","tid","pos","density")

    # add features info
    dt <- dt |>
      dplyr::filter(density > 0) |>
      dplyr::filter(gene_name %in% select_gene$gene_name) |>
      dplyr::left_join(y = ganno,by = "tid",multiple = "all") |>
      dplyr::filter(pos >= utr5 + upstreamNt & pos <= utr5 + cds - downstreamNt)

    # calculate total density for gene
    total_density <- dt |>
      dplyr::group_by(tid) |>
      dplyr::summarise(sum_density = sum(density))

    # calculate wi and pi for each gene
    if(is.null(group)){
      group <- sam_file
    }

    ps_df <- dt |>
      dplyr::left_join(y = total_density,by = "tid",multiple = "all") |>
      dplyr::mutate(abs_pos = pos - utr5,
                    wi = (2*abs_pos -(cds + 1))/(cds - 1),
                    pi = density*wi/sum_density) |>
      dplyr::group_by(gene_name) |>
      dplyr::summarise(sum_pi = sum(pi)) |>
      dplyr::mutate(group = group,sample = sam_file)
  }

  # ========================================================================================
  # loop to calclate polarity scores for each samples
  # ========================================================================================
  # x = 1
  purrr::map_df(seq_along(input_file),function(x){
    calc_PS_fun(gtf_file = gtf_file,
                input_file = input_file[x],
                density_file = density_file[x],
                minCounts = minCounts,
                upstreamNt = upstreamNt,
                downstreamNt = downstreamNt,
                group = group[x])
  }) -> ps_df

  return(ps_df)
}





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
#'   Default is "grey50".
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
                                orf_col = "grey50",
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
