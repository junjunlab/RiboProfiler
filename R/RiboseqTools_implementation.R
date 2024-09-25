globalVariables(c("gene", "hi_CI", "lo_CI", "mean_alpha", "median", "next.hi", "next.lo", "overlap",
                  "overlap.xmax", "overlap.xmin", "overlap.ymax", "overlap.ymin", "setDT", "winmid", "xmax",
                  "xmin","codon_pos", "density_sd", "ratio", "boot"))

.defaults <- getFromNamespace(".defaults","RiboSeqTools")
calc_total_counts <- getFromNamespace("calc_total_counts","RiboSeqTools")
set_excluded <- getFromNamespace("set_excluded","RiboSeqTools")
set_total_counts <- getFromNamespace("set_total_counts","RiboSeqTools")
set_downsampled <- getFromNamespace("set_downsampled","RiboSeqTools")
get_default_param <- getFromNamespace("get_default_param","RiboSeqTools")




trans_experiment <- function(..., object = NULL,
                             .bin = c('bynuc', 'byaa'), .exclude = NULL) {
  bin <- match.arg(.bin, several.ok=TRUE)
  ret <- mapply(function(...) {
    paths <- rlang::list2(...)
    ret <- sapply(paths, function(path) {
      ret <- list()

      m <- transform_to_sparse_matrix(object = object,
                                      sample_name = path)

      libsize <- sum(m, na.rm=TRUE)
      if ('bynuc' %in% .bin)
        ret$bynuc <- m
      if ('byaa' %in% .bin) {
        mod <- ncol(m) %% 3
        if (mod) {
          rlang::warn(sprintf("Length of at least one CDS in %s is not divisible by 3", path))
          m <- m[,1:(ncol(m) - mod)]
        }
        ret$byaa <- m[,seq(1, ncol(m), by=3), drop=FALSE] + m[,seq(2, ncol(m), by=3), drop=FALSE] + m[,seq(3, ncol(m), by=3), drop=FALSE]
      }
      if (!is.null(.exclude) && length(.exclude) > 0)
        ret <- lapply(ret, function(x) x[!(rownames(x) %in% .exclude),])
      ret
    }, simplify=FALSE)
    ret
  },  ..., SIMPLIFY=FALSE)
  names(ret) <- 1:unique(purrr::map_int(list(...), length))
  ret
}






#' Transform ribosomeObj Object of shifted.data into sparsed matrix
#'
#' @param object ribosomeObj Object.
#' @param sample_name sample name.
#'
#' @return A sparsed matrix containing gene read counts on cds region.
#' @export
transform_to_sparse_matrix <- function(object = NULL,sample_name = NULL){
  shifted_counts <- object@shifted.data

  data.table::setDT(shifted_counts)

  ref <- object@longest.annotation
  ref <- ref[,c("gene_name","tid","utr5","cds")]

  tmp <- shifted_counts[sample == sample_name, .(counts = sum(counts)), by = .(trans_id, trans_pos)]

  tmp <- tmp %>%
    dplyr::left_join(y = ref,by = c("trans_id" = "tid")) %>%
    dplyr::mutate(trans_pos = trans_pos - utr5) %>%
    dplyr::filter(trans_pos > 0 & trans_pos <= cds & gene_name != "noGeneName") %>%
    dplyr::select(-utr5,-trans_id)

  gene_info <- tmp[,c("gene_name","cds")] %>% unique()

  gene_info_all_len <- purrr::map2_dfr(gene_info$gene_name,gene_info$cds,
                                       function(gene_name,cds_len){
                                         data.frame(gene_name = gene_name,trans_pos = 1:cds_len)
                                       })


  # all position counts
  all_pos_counts <- gene_info_all_len %>%
    dplyr::left_join(y = tmp[,-c("cds")],
                     by = c("gene_name","trans_pos"),
                     relationship = "many-to-many") %>%
    tidyr::replace_na(list(counts = 0))

  # longer to wider
  df_wide <- data.table::dcast(setDT(all_pos_counts),
                               gene_name ~ trans_pos,
                               value.var = "counts",
                               drop = FALSE,fill = 0) %>%
    tibble::column_to_rownames(var = "gene_name")

  acols <- ncol(df_wide) %% 3
  if (acols > 0) {
    acols <- ncol(df_wide):(ncol(df_wide) - acols + 1)
    if (all(is.na(df_wide[,acols])))
      df_wide <- df_wide[,-acols]
  }

  # turn into sparsed matrix
  df_wide_sparseed <- Matrix::Matrix(data = as.matrix(df_wide), sparse = TRUE)

  return(df_wide_sparseed)
}






#' load selective ribosome data
#'
#' @param ... Name-value pairs of lists. The name of each argument will be the name of an experiment.
#'      The name of each element will be the sample type (e.g. TT for total translatome), the value
#'      of each element must be a character vector of file paths, where each file is a read count
#'      table of one replicate experiment. Replicate order must match between sample types.
#' @param object object ribosomeObj Object.
#' @param normalize Normalize the read counts to library size? Output will then be in RPM.
#' @param bin Bin the data. \code{bynuc}: No binning (i.e. counts per nucleotide). \code{byaa}: Bin by residue.
#' @param exclude Genes to exclude in all future analyses. This genes will also be excluded from total read count
#'      calculation. Note that the raw count tables will not be modified. Named list with names corresponding to
#'      experiments. If a character vector of gene names is given, these genes will be excluded from all
#'      experiments.
#' @param defaults Default parameters of the data set.
#' @return An object of class \code{serp_data}
#'
#' @importFrom stats median
#'
#' @export
load_selective_ribosome <- function(..., object = NULL,
                                    normalize = FALSE,
                                    bin = c('bynuc', 'byaa'),
                                    exclude = list(), defaults = list()) {
  experiments <- rlang::list2(...)
  if (is.null(names(experiments)) || sum(nchar(names(experiments)) == 0) > 0){
    rlang::abort("all experiments must be named")
  }


  bin <- match.arg(bin,c('bynuc', 'byaa'))

  if (!is.list(exclude) && is.character(exclude)){
    exclude <- purrr::map(experiments, function(...)exclude)
  }

  data <- sapply(experiments, purrr::lift_dl(trans_experiment),
                 object = object, .bin=bin, simplify=FALSE)

  what <- ifelse('bynuc' %in% bin, 'bynuc', 'byaa')

  # load gene annotation
  ref <- object@longest.annotation
  ref <- ref[,c("gene_name","cds")]
  colnames(ref) <- c("gene","length")

  ref$cds_length <- ref$length %/% 3

  ret <- list(ref = ref, data = data, normalized = FALSE)
  ret <- structure(ret, class="serp_data")
  ret <- set_excluded(ret, exclude)
  ret <- set_total_counts(ret, calc_total_counts(ret, what))

  if (normalize){
    ret <- normalize(ret)
  }

  defaults <- purrr::list_modify(.defaults, !!!defaults)

  if (defaults$bin == 'byaa' && !('byaa' %in% bin))
    defaults$bin <- 'bynuc'

  ret$defaults <- defaults

  set_downsampled(ret, FALSE)
}







#' Plot a profile along a gene
#'
#' For each position, a confidence interval is plotted. Transparency reflects the total number of reads
#' contributing to the confidence interval.
#'
#' @param object object ribosomeObj Object.
#' If \code{type} is \code{enrichment}, a binomial confidence interval for enrichment of \code{sample1}
#' compared to \code{sample2} is calculated using \code{\link{binom_ci_profile}}. The \code{samples}
#' argument is ignored.
#'
#' If \code{type} is \code{samples}, a Poisson confidence interval for RPM in samples \code{samples}
#' is calculated using \code{\link{pois_ci_profile}}. The upper and lower bounds are then divided by
#' \code{window_size} to indicate confidence in the local smoothed read density. The \code{sample1}
#' and \code{sample2} arguments are ignored.
#'
#' @param object A \code{serp_data} object. Must contain raw (unnormalized) read counts.
#' @param gene_name Name of the gene/ORF to plot.
#' @param type Plot type. One of \code{enrichment}, \code{rpm}.
#' @param samples Samples to plot. If missing, all samples will be plotted.
#' @param sample1 Name of the first sample (the numerator). If missing, the default sample1 of the data set
#'      will be used.
#' @param sample2 Name of the second sample (the denominator). If missing, the default sample2 of the data set
#'      will be used.
#' @param exp Character vector of experiments to plot. If missing, all experiments are plotted.
#' @param rep Character vector of replicates to plot. If missing, all replicates will be plotted.
#' @param bin Bin level (\code{bynuc} or \code{byaa}). If missing, the default binning level of the data set
#'      will be used.
#' @param window_size Neighborhood size for the confidence interval calculation in nucleotides. If missing, the default
#'      window size of the data set will be used.
#' @param conf.level Confidence level.
#'
#' @return A data.frame.
#' @export
get_serp_plot_df <- function(object = NULL,
                             type = c("enrichment","rpm"),
                             gene_name = NULL,
                             samples = NULL,
                             sample1 = NULL, sample2 = NULL,
                             exp, rep,
                             bin = c("byaa","bynuc"),
                             window_size = 45, conf.level = 0.95){
  bin <- match.arg(bin,c("byaa","bynuc"))
  type <- match.arg(type,c("enrichment","rpm"))
  cds_len <- subset(object$ref,gene == gene_name)$cds_length

  # check type
  if (type == "enrichment") {
    sample1 <- get_default_param(data, sample1)
    sample2 <- get_default_param(data, sample2)

    df <- RiboSeqTools::binom_ci_profile(data = object, gene = gene_name,
                                         sample1 = sample1, sample2 = sample2,
                                         exp, rep ,
                                         bin = bin,
                                         window_size = window_size, conf.level=conf.level) %>%
      dplyr::mutate(alpha=1/((1/!!rlang::sym(paste0('win_', sample1)) + 1/!!rlang::sym(paste0('win_', sample2))) * window_size))

    dfgroups <- list(rlang::sym("exp"), rlang::sym("rep"))
  } else if (type == "rpm") {
    df <- RiboSeqTools::pois_ci_profile(data = object, gene = gene_name,
                                        samples = samples,
                                        exp , rep ,
                                        bin = bin,
                                        window_size = window_size, conf.level=conf.level) %>%
      # dplyr::mutate(alpha=counts / window_size, lo_CI=lo_CI / window_size, hi_CI=hi_CI / window_size) +
      dplyr::mutate(alpha=counts / window_size)

    dfgroups <- list(rlang::sym("exp"), rlang::sym("rep"), rlang::sym("sample"))
  }

  df_plot <- dplyr::group_by(df, !!!dfgroups) %>%
    dplyr::mutate(xmin = (winmid - 0.5) ,
                  xmax = (winmid + 0.5),
                  next.lo = dplyr::lead(lo_CI, default=0),
                  next.hi = dplyr::lead(hi_CI, default=0),
                  overlap = dplyr::if_else(lo_CI - next.hi > 0, 1L, dplyr::if_else(next.lo - hi_CI > 0, 2L, 0L)),
                  overlap.xmin = xmax - (xmax - xmin) * 0.5,
                  overlap.xmax = dplyr::lead(xmin) + (dplyr::lead(xmax) - dplyr::lead(xmin)) * 0.5,
                  overlap.ymin = dplyr::recode(overlap, `0`=NA_real_, `1`=next.hi, `2`=hi_CI),
                  overlap.ymax = dplyr::recode(overlap, `0`=NA_real_, `1`=lo_CI, `2`=next.lo),
                  mean_alpha = mean(c(alpha, dplyr::lead(alpha)), na.rm=TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene_name = gene_name)

  attr(df_plot,"bin") <- bin
  attr(df_plot,"type") <- type
  attr(df_plot,"cds_len") <- cds_len

  return(df_plot)
}








#' Single-gene Ribosome Profiling Plot
#'
#' This function creates a visualization for ribosome profiling of a single gene,
#' supporting clustering, functional annotation, and publication-ready plots.
#' The input data can be transformed and structured into a faceted plot,
#' with optional inclusion of coding sequence (CDS) annotations.
#'
#' @param data A dataframe containing ribosome profiling data. The attribute "bin"
#'   should indicate the bin type (`"byaa"` or `"bynuc"`), and the attribute "type"
#'   should determine if the data represents `"enrichment"` or ribosome density.
#' @param transform Character string determining the scale transformation for the y-axis.
#'   Default is `"log2"`.
#' @param aes_color Aesthetics mapping for color, usually representing a grouping variable
#'   like `"sample"` or `"rep"`. Default is `"sample"`.
#' @param y_threshold Numeric threshold for drawing a horizontal dashed line when the
#'   plot type is enrichment. Default is `1`.
#' @param shadow_alpha_range A numeric vector of length 2 controlling the range of the
#'   alpha transparency for shadow regions. Default is `c(0.2,1)`.
#' @param facet_grid_params A list of additional parameters passed to `facet_grid`.
#'   Default is an empty list.
#' @param add_cds_structure Logical, whether to add a coding sequence (CDS) annotation
#'   structure to the plot. Default is `TRUE`.
#' @param cds_col Color of the CDS annotation line. Default is `"grey"`.
#' @param cds_anno_col Color of the CDS annotation labels. Default is `"#996600"`.
#' @param cds_region_width Numeric, the width of the CDS annotation line. Default is `5`.
#' @param cds_region_scale_x Numeric, scaling factor for the x-axis of the CDS region
#'   annotation. Default is `NULL` (no scaling).
#'
#' @return A ggplot object with ribosome profiling visualizations, optionally including
#'   structure annotations like the CDS.
#'
#'
#' @import ggplot2 dplyr ggside
#' @export
single_gene_profile <- function(data = NULL,
                                transform = "log2",
                                aes_color = "sample",
                                y_threshold = 1,
                                shadow_alpha_range = c(0.2,1),
                                facet_grid_params = list(),
                                add_cds_structure = TRUE,
                                cds_col = "grey",
                                cds_anno_col = "#996600",
                                cds_region_width = 5,
                                cds_region_scale_x = NULL){
  # add group
  data <- data %>%
    dplyr::mutate(group = paste(rep,sample,sep = "-"))

  # check bin type
  if(attr(data,"bin") == "byaa"){
    xlabel <- "Ribosome along transcript position\n (codons / amino acids)"
  }else{
    xlabel <- "Ribosome along transcript position\n (nucleotides)"
  }

  # check plot type
  if(attr(data,"type") == "enrichment"){
    ylabel <- "Enrichment (IP / Total)"
    y <- "ratio_mean"
    hline <- geom_hline(yintercept = y_threshold,lty = "dashed",color = "black")
  }else{
    ylabel <- "Ribosome density (RPM)"
    y <- "RPM_mean"
    hline <- NULL
  }


  color <- rlang::ensym(aes_color)
  y <- rlang::ensym(y)

  # ============================================================================
  # plot
  # ============================================================================
  pmain <-
    ggplot(data = data) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = lo_CI, ymax = hi_CI,
                  alpha = alpha,
                  fill = !!color),color = NA,linewidth = 0,show.legend = F) +
    geom_rect(aes(xmin = overlap.xmin, xmax = overlap.xmax,
                  ymin = overlap.ymin, ymax = overlap.ymax,
                  alpha = mean_alpha,
                  fill = !!color),
              show.legend = F,color = NA,linewidth = 0,
              data = function(y) dplyr::filter(y, overlap > 0, !is.na(overlap.xmax))) +
    geom_line(aes(x = winmid,y = !!y,color = !!color,group = group)) +
    hline +
    theme_bw() +
    # facet_grid(cols = vars(gene_name),rows = vars(rep)) +
    do.call(facet_grid,modifyList(list(cols = vars(gene_name),
                                       rows = vars(rep)),
                                  facet_grid_params)) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.placement = "outside",
          strip.text = element_text(face = "bold",size = rel(1)),
          strip.text.x = element_text(face = "bold.italic",size = rel(1)),
          strip.background = element_blank(),
          strip.text.y.right = element_text(angle = 0,hjust = 0)) +
    scale_alpha_continuous(transform ="sqrt",
                           # breaks = c(1,3,10,30),
                           limits = c(0,30),
                           name = "prec", range = shadow_alpha_range) +
    scale_y_continuous(transform = transform) +
    ylab(ylabel) +
    xlab(xlabel)

  # ==========================================================================
  # add structure
  # ==========================================================================
  if(add_cds_structure == TRUE){
    cds_rg_struc <- data.frame(x = 0,xend = attr(data,"cds_len"),
                               y = 1,yend = 1)

    # plot
    pmain_anno <- pmain +
      ggside::geom_xsidesegment(data = cds_rg_struc,
                                mapping = aes(x = x,xend = xend,y = y,yend = yend),
                                linewidth = cds_region_width,color = cds_col) +
      theme(ggside.panel.background = element_blank(),
            ggside.panel.border = element_blank(),
            ggside.panel.scale.x = cds_region_scale_x) +
      ggside::scale_xsidey_continuous(breaks = NULL) +
      ggside::ggside(x.pos = "top",collapse = "x")
  }else{
    pmain_anno <- pmain
  }

  return(pmain_anno)
}







#' Generate Metagene Profile Plot
#'
#' This function creates a metagene profile plot using `ggplot2`, visualizing
#' ribosome profiling data or other biological signals aligned to a specific reference
#' point such as the start or stop codon. It also supports confidence intervals (CI)
#' for bootstraps and faceting for additional dimensions of group or experiment.
#'
#' @param df A dataframe containing positional metagene data. Should include columns
#'   like `pos` (position), `summary` (signal summary), and a grouping variable for fill/color.
#' @param ytrans A string specifying the transformation for the y-axis. Default is `"log2"`.
#'   Example values are `"log10"`, `"sqrt"` or `"identity"` for no transformation.
#' @param align A string indicating how the metagene data is aligned. Can be `"start"`
#'   for start codon or `"stop"` for stop codon. Used in formatting the x-axis label. Default is `"start"`.
#' @param aes_col Aesthetic mapping for the fill and color of the plot. This variable
#'   represents the grouping factor for different metagene profiles. Should match a column name in `df`.
#' @param facet_grid_params A list of parameters to pass to `facet_grid()`, allowing for faceted plots by different factors (e.g., experimental replicates).
#'   Default is an empty list.
#' @param conf.level The confidence level for the bootstrap confidence intervals (CI). Must be a number between `0` and `1`. Default is `0.95`.
#' @param ci.alpha The alpha transparency value (between `0` and `1`) for the confidence interval shade. Default is `0.3`.
#'
#' @return A ggplot object visualizing the metagene profile with confidence intervals.
#'
#' @importFrom stats quantile
#'
#' @export
metagene_profile <- function(df = NULL,
                             ytrans = "log2",
                             align = "start",
                             aes_col = NULL,
                             facet_grid_params = list(),
                             conf.level = 0.95, ci.alpha = 0.3){
  col <- rlang::ensym(aes_col)

  # facet
  if(length(facet_grid_params) > 0){
    facet_layer <- do.call(facet_grid,modifyList(list(),facet_grid_params))
  }else{
    facet_layer <- NULL
  }

  # plot
  ggplot(data = df, aes(x = pos, y = summary,fill = !!col,color = !!col)) +
    stat_summary(data = function(x) dplyr::filter(x, boot),
                 aes(color = NULL),
                 geom = 'ribbon',
                 fun.ymin = function(x) quantile(x, 0.5 * (1 - conf.level)),
                 fun.ymax = function(x) quantile(x, 1 - 0.5 * (1 - conf.level)),
                 alpha = ci.alpha) +
    geom_line(aes(y = summary), data = function(x) dplyr::filter(x, !boot)) +
    scale_y_continuous(transform = ytrans) +
    labs(x = sprintf('distance from %s codons', align), y = ylab) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)), color = FALSE) +
    theme_bw() +
    facet_layer +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"),
          strip.placement = "outside",
          strip.text = element_text(face = "bold",size = rel(1)),
          strip.text.x = element_text(face = "bold.italic",size = rel(1)),
          strip.background = element_blank(),
          strip.text.y.right = element_text(angle = 0,hjust = 0))
}
