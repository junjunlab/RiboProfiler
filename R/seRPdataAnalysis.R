globalVariables(c("codon_pos","density_sd","mean_ratio", "nt_pos", "ratio"))

#' pre_enrichment_data function
#'
#' This function performs enrichment analysis on ribosome profiling data.
#'
#' @param riboIP_file A character vector containing the file path(s) to the
#' ribosome profiling IP sample(s).
#' @param riboInput_file A character vector containing the file path(s) to the
#' ribosome profiling input control sample(s).
#' @param output_file A character vector containing the desired output file
#' path(s) for each sample's enrichment analysis results.
#' @param pythonPath An optional parameter specifying the file path to a specific
#' Python installation to use (if not already set as the default in reticulate).
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' pre_enrichment_data(riboIP_file = "path/to/IP_sample.bam",
#'                     riboInput_file = "path/to/Input_sample.bam",
#'                     output_file = "path/to/output_file.txt")
#' }
#'
#' @export
pre_enrichment_data <- function(riboIP_file = NULL,
                                riboInput_file = NULL,
                                output_file = NULL,
                                pythonPath = NULL){
  # check python
  reticulate::py_config()
  if(reticulate::py_available() == FALSE){
    message("Please install python first!")
  }else{
    if(!is.null(pythonPath)){
      reticulate::use_python(pythonPath)
    }

    # run code
    pyscript.path = system.file("extdata", "enrichmentAnalysis.py", package = "RiboProfiler")
    reticulate::source_python(pyscript.path)

    # loop
    lapply(seq_along(riboIP_file), function(x){
      reticulate::py$enrichmentAnalysis(IPfile = riboIP_file[x],
                                        Inputfile = riboInput_file[x],
                                        outputFile = output_file[x])

      # print
      message(paste(output_file[x]," has been processed!",sep = ""))
    }) -> tmp
  }
}


#' Preprocess ribosome profiling data using sliding window enrichment analysis.
#'
#' This function preprocesses ribosome profiling data by running sliding window
#' enrichment analysis (SWEA) on the input files.
#'
#' @param riboIP_file A character vector of file paths to the ribosome profiling
#' IP samples.
#' @param riboInput_file A character vector of file paths to the ribosome profiling
#' Input samples.
#' @param gene_list A character vector of file path to the gene list.
#' @param output_file A character vector of output file names for SWEA results.
#' @param pythonPath A character string specifying the path to the Python executable.
#' Defaults to NULL, which results in automatic detection of an available Python
#' installation.
#'
#' @details This function calls a Python script provided with the RiboProfiler
#' package to perform sliding window enrichment analysis on ribosome profiling data.
#'
#' @return The function invisibly returns a list of processed files.
#'
#' @examples
#' \dontrun{
#' # Preprocess ribosome profiling data using sliding window enrichment analysis
#' pre_slideWindow_enrichment_data(riboIP_file = c("IP1.density.txt", "IP2.density.txt"),
#'                                 riboInput_file = c("Input1.density.txt", "Input2.density.txt"),
#'                                 gene_list = c("gene1","gene2"),
#'                                 output_file = c("enrich_IP1.txt", "enrich_IP2.txt"))
#' }
#'
#' @export
pre_slideWindow_enrichment_data <- function(riboIP_file = NULL,
                                            riboInput_file = NULL,
                                            gene_list = NULL,
                                            output_file = NULL,
                                            pythonPath = NULL){
  # check python
  reticulate::py_config()
  if(reticulate::py_available() == FALSE){
    message("Please install python first!")
  }else{
    if(!is.null(pythonPath)){
      reticulate::use_python(pythonPath)
    }

    # run code
    pyscript.path = system.file("extdata", "slidewindowEnrichmentAnalysis.py", package = "RiboProfiler")
    reticulate::source_python(pyscript.path)

    # loop
    lapply(seq_along(riboIP_file), function(x){
      reticulate::py$batchEnrichment(inputIPFile = riboIP_file[x],
                                     inputInputFile = riboInput_file[x],
                                     geneList = gene_list,
                                     ouputFile = output_file[x])

      # print
      message(paste(output_file[x]," has been processed!",sep = ""))
    }) -> tmp
  }
}


#' Load and enrich ribosome profiling data
#'
#' This function processes raw or smoothed ribosome profiling data and returns
#' a data frame with mean density and standard deviation values for each codon
#' position, gene, sample, and type (ribosome profiling).
#'
#' @param data_type The type of input data: "raw" or "sw".
#' @param enrich_data A character vector of file paths to the input data files.
#' @param gene_names A character vector of gene names to filter on.
#' @param sample_name A character vector of sample names to use instead of file names.
#'
#' @return A data frame with columns for gene_name, trans_id, codon_pos, type,
#' sample, density, and (density_sd).
#'
#' @examples
#' \dontrun{
#' load_enrich_data(data_type = "raw", enrich_data = c("sample1.txt", "sample2.txt"),
#'                  gene_names = c("GAPDH", "ACTB"))
#' load_enrich_data(data_type = "sw", enrich_data = c("sample1.txt", "sample2.txt"),
#'                  sample_name = c("Sample A", "Sample B"))
#' }
#'
#' @export
load_enrich_data <- function(data_type = c("raw","sw"),
                             enrich_data = NULL,
                             gene_names = NULL,
                             sample_name = NULL){
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  data_type <- match.arg(data_type,c("raw","sw"))
  # ============================================================================
  # 1_raw density
  # ============================================================================
  # process data
  if(data_type == "raw"){
    # loop read
    # x = 1
    plyr::ldply(1:length(enrich_data), function(x){
      tmp <- vroom::vroom(enrich_data[x],
                          delim = "\t",show_col_types = FALSE,
                          col_names = c('id','pos','ip_denisty','input_density','ratio'),
                          col_select = c("id","pos","ratio")) |>
        data.table::as.data.table()

      # add gene name and filter target gene
      tmp <- tmp[, c("gene_name","trans_id") := data.table::tstrsplit(id, "|", fixed = TRUE)[c(1,3)]] |>
        dplyr::filter(gene_name %in% gene_names)

      # make trans pos into codon pos
      # x = "ARC1"
      plyr::ldply(unique(tmp$gene_name),function(x){
        tmp1 <- tmp[gene_name == x]
        start <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',4) %>% as.numeric()
        end <- sapply(strsplit(tmp1$id[1],split = '\\|'),'[',5) %>% as.numeric()

        # 1.add nt pos
        tmp1$nt_pos <- tmp1$pos - start + 1

        # 2.transform to codon pos
        sq <- seq(1,(end - start + 1),3);rg <- c(1:length(sq))
        plyr::ldply(1:length(sq),function(z){
          tmp2 = tmp1[nt_pos >= sq[z] & nt_pos <= sq[z] + 2
          ][,.(mean_ratio = mean(ratio)),by = .(id,gene_name,trans_id)
          ][,`:=`(codon_pos = rg[z])]
          return(tmp2)
        }) -> codon_res

        # 3.add to continues codon positions
        # s = 375
        plyr::ldply(rg,function(s){
          if(s %in% codon_res$codon_pos){
            tmp3 <- codon_res[which(codon_res$codon_pos == s),]
          }else{
            tmp3 <- data.frame(id = codon_res$id[1],
                               gene_name = codon_res$gene_name[1],
                               trans_id = codon_res$trans_id[1],
                               mean_ratio = 0,codon_pos = s)
          }
          return(tmp3)
        }) -> continues_codon_res
        return(continues_codon_res)
      }) -> final_res

      # add sample info
      if(is.null(sample_name)){
        sample <- enrich_data[x]
      }else{
        sample <- sample_name[x]
      }
      final_res$type <- "ribo"
      final_res$sample <- sample

      return(final_res)
    }) -> df_ratio

    # mean for replicates
    merge_rep <- df_ratio |>
      dplyr::group_by(gene_name,trans_id,codon_pos,type,sample) |>
      dplyr::summarise(mean_rep_ratio = mean(mean_ratio),
                       mean_sd = stats::sd(mean_ratio))
    merge_rep[is.na(merge_rep)] <- 0
    colnames(merge_rep)[6:7] <- c("density","density_sd")
  }else if(data_type == "sw"){
    # ============================================================================
    # 2_window smooth
    # ============================================================================
    # loop read
    # x = 1
    plyr::ldply(1:length(enrich_data), function(x){
      tmp <- vroom::vroom(enrich_data[x],
                          delim = "\t",show_col_types = FALSE,
                          col_names = c('id','codon_pos','mean_ratio')) |>
        data.table::as.data.table()

      # add gene name and filter target gene
      tmp <- tmp[, c("gene_name","trans_id") := data.table::tstrsplit(id, "|", fixed = TRUE)[c(1,3)]]

      # filter gene
      if(!is.null(gene_names)){
        final_res <- tmp %>%
          dplyr::filter(gene_name %in% gene_names)
      }else{
        final_res <- tmp
      }

      # add sample info
      if(is.null(sample_name)){
        sample <- enrich_data[x]
      }else{
        sample <- sample_name[x]
      }

      final_res$type <- "ribo"
      final_res$sample <- sample

      return(final_res)
    }) -> df_ratio

    # mean for replicates
    merge_rep <- df_ratio |>
      dplyr::group_by(gene_name,trans_id,codon_pos,type,sample) |>
      dplyr::summarise(mean_rep_ratio = mean(mean_ratio),
                       mean_sd = stats::sd(mean_ratio))
    merge_rep[is.na(merge_rep)] <- 0
    colnames(merge_rep)[6:7] <- c("density","density_sd")
  }

  return(merge_rep)
}
