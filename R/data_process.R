globalVariables(c("average_exp", "norm_counts", "trancript_len"))

#' Get information about longest transcripts in a GTF file
#'
#' This function uses a Python script to extract information about the longest
#' transcripts in a GTF file. The script generates an output file containing the
#' transcript ID, gene ID,... and length of the longest transcript for each gene.
#'
#' @param gtf_file path to the GTF file, ".gz" compressed format is acceptted.
#' @param out_file path to the output file.
#' @param pythonPath path to the Python executable (optional).
#'
#' @return This function does not return anything. It writes the output to a file.
#'
#' @importFrom reticulate py_config py_available use_python py_module_available
#' py_install source_python
#'
#' @examples
#' \dontrun{
#' pre_longest_trans_info(gtf_file = "my_gtf_file.gtf",
#'                        out_file = "longest_transcripts.txt",
#'                        pythonPath = "/usr/bin/python3")
#'}
#'
#' @export
pre_longest_trans_info <- function(gtf_file = NULL,
                                   out_file = NULL,
                                   pythonPath = NULL){
  reticulate::py_config()
  if(reticulate::py_available() == FALSE){
    message("Please install python first!")
  }else{
    if(!is.null(pythonPath)){
      reticulate::use_python(pythonPath)
    }

    # check pyfaidx
    if (!reticulate::py_module_available("gzip")) {
      cat("Installing gzip ...\n")
      reticulate::py_install("gzip")
    }

    # import gzip
    tryCatch({
      py <- reticulate::import("gzip")
      # print(py)
    }, error = function(e) {
      cat("Error: gzip is not available.\n")
    })

    # sort gtf first
    gtf <- rtracklayer::import.gff(gtf_file,format = "gtf")

    # sort by strand
    gtf_plus <- gtf |> data.frame() |>
      dplyr::filter(strand %in% "+") |>
      dplyr::arrange(seqnames,gene_name,gene_id,transcript_id,type,start,end)

    gtf_neg <- gtf |> data.frame() |>
      dplyr::filter(strand %in% "-") |>
      dplyr::arrange(seqnames,gene_name,gene_id,transcript_id,type,dplyr::desc(start),dplyr::desc(end))

    sorted_gtf <- rbind(gtf_plus,gtf_neg)

    # output
    output_name = paste(gtf_file,".sorted.gtf",sep = "")
    rtracklayer::export.gff(sorted_gtf,
                            con = output_name,
                            format = "gtf")

    # run code
    pyscript.path = system.file("extdata", "getLongestTransInfo.py", package = "RiboProfiler")
    reticulate::source_python(pyscript.path)
    suppressMessages(
      reticulate::py$getLongestTransInfo(gtf_file = output_name,
                                         longest_file = out_file)
    )
  }
}







#' Calculate RNA coverage data
#'
#' This function generates pre-processing RNA coverage data for a given SAM/BAM file.
#'
#' @param sam_file A character vector specifying the path(s) of the input SAM file(s).
#' @param bam_file If using sam file runs slowly, try using bam file instead.
#' @param out_file A character vector specifying the name(s) of the output file(s).
#' The output file(s) will be saved in the "2.density-data" directory.
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Generate pre-processing RNA coverage data for a single SAM file
#' pre_rna_coverage_data(sam_file = "/path/to/sam_file.sam", out_file = "coverage_data.txt")
#'
#' # Generate pre-processing RNA coverage data for multiple SAM files
#' pre_rna_coverage_data(sam_file = c("/path/to/sam_file_1.sam", "/path/to/sam_file_2.sam"),
#'                        out_file = c("coverage_data_1.txt", "coverage_data_2.txt"))
#' }
#'
#' @export
pre_rna_coverage_data <- function(sam_file = NULL,
                                  bam_file = NULL,
                                  out_file = NULL){
  if(!dir.exists("2.density-data")){
    dir.create("2.density-data")
  }

  # check input file type
  if(!is.null(bam_file)){
    lapply(seq_along(bam_file), function(x){
      outFile_tmp = paste("2.density-data/",out_file[x],sep = "")

      # bam total mapped reads
      tryCatch(expr = Rsamtools::idxstatsBam(bam_file[x]),
               error = function(e){
                 Rsamtools::indexBam(files = bam_file)
               })

      total_mapped_reads <- sum(Rsamtools::idxstatsBam(bam_file[x])$mapped)

      # params
      scanbamparam <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
      pileupparam <- Rsamtools::PileupParam(max_depth = 10000000,
                                            distinguish_strands = FALSE,
                                            distinguish_nucleotides = FALSE)

      # get depth
      depth <- Rsamtools::pileup(file = bam_file[x],
                                 pileupParam = pileupparam,
                                 scanBamParam = scanbamparam )

      # rpm normalization
      depth$rpm <- (depth$count / total_mapped_reads)*10^6

      # output
      data.table::fwrite(depth,file = outFile_tmp,sep = "\t",col.names = FALSE,nThread = parallel::detectCores())

      message(paste(bam_file[x]," has been processed!",sep = ""))
    }) -> tmp
    return(NULL)
  }else if(!is.null(sam_file)){
    JuliaCall::julia_setup(installJulia = TRUE)

    JuliaCall::julia_library("XAM")

    script_path <- paste0('include("',
                          system.file("extdata", "CalculateRNACoverage.jl",
                                      package = "RiboProfiler"),
                          '")',collapse = "")

    calculateRNACoverage <- JuliaCall::julia_eval(script_path)

    # excute function
    lapply(seq_along(sam_file), function(x){
      outFile_tmp = paste("2.density-data/",out_file[x],sep = "")
      calculateRNACoverage(inputFile = sam_file[x],
                           outputFile = outFile_tmp,
                           type = "coverage")
      message(paste(sam_file[x]," has been processed!",sep = ""))
    }) -> tmp
    return(NULL)
  }
}


#' Calculate density data for a given gene annotation file
#'
#' This function calculates density data for a given gene annotation file using
#' the XAM package in Julia. It writes the output files to the "2.density-data"
#' directory.
#'
#' @param gene_anno A gene annotation file in the format accepted by XAM
#' (default is NULL).
#' @param density_file A character vector of input file names for which to
#' calculate density data (default is NULL).
#' @param out_file A character vector of output file names corresponding to each
#' input file (default is NULL). The length of this vector must be equal to the
#' length of \code{density_file}.
#'
#' @return This function does not return anything; it writes the output files to
#' disk.
#'
#' @importFrom JuliaCall julia_eval julia_library julia_setup
#'
#' @examples
#' \dontrun{
#' # Calculate density data for a single input file
#' pre_gene_trans_density(gene_anno = "gene_annotation.txt",
#'                        density_file = "input_file.txt",
#'                        out_file = "output_file.txt")
#'
#' # Calculate density data for multiple input files
#' pre_gene_trans_density(gene_anno = "gene_annotation.txt",
#'                        density_file = c("input_file1.txt", "input_file2.txt"),
#'                        out_file = c("output_file1.txt", "output_file2.txt"))
#' }
#'
#' @export
pre_gene_trans_density <- function(gene_anno = NULL,
                                   density_file = NULL,
                                   out_file = NULL){
  if(!dir.exists("2.density-data")){
    dir.create("2.density-data")
  }

  JuliaCall::julia_setup(installJulia = TRUE)

  script_path <- paste0('include("',
                        system.file("extdata", "GetGeneSinglePosDensity.jl",
                                    package = "RiboProfiler"),
                        '")',collapse = "")

  getGeneSinglePosDensity <- JuliaCall::julia_eval(script_path)

  # excute function
  lapply(seq_along(density_file), function(x){
    inputFile_tmp = paste("2.density-data/",density_file[x],sep = "")
    outFile_tmp = paste("2.density-data/",out_file[x],sep = "")

    getGeneSinglePosDensity(geneInfo = gene_anno,
                            inputFile = inputFile_tmp,
                            outputFile = outFile_tmp)
    message(paste(inputFile_tmp," has been processed!"))
  }) -> tmp
  return(NULL)
}






#' Load gene expression data from ribo and rna files
#'
#' This function reads in gene expression data from ribo and rna files, filters
#' by a given gene list, and merges them into a single data frame. The resulting
#' data frame includes columns for gene name, transcript ID, position, density,
#' sample name, and type of expression data (ribo or rna).
#'
#' @param mapping_type The mapping type for your sam files, "genome" or "transcriptome".
#' @param ribo_file A character vector specifying the names of ribo file(s).
#' @param rna_file A character vector specifying the names of rna file(s).
#' @param sample_name A character vector specifying the name(s) of the sample(s).
#' @param gene_list A character vector specifying the gene names to be included.
#'
#' @return A data frame containing gene expression data from both ribo and rna
#' files, filtered by the specified gene list and merged into a single data frame.
#'
#' @examples
#' \dontrun{
#' load_track_data(ribo_file = c("file1_ribo.txt", "file2_ribo.txt"),
#'                 rna_file = c("file1_rna.txt", "file2_rna.txt"),
#'                 sample_name = c("sample1", "sample2","sample1", "sample2"),
#'                 gene_list = c("gene1", "gene2"))
#' }
#'
#' @export
load_track_data <- function(mapping_type = c("genome","transcriptome"),
                            ribo_file = NULL,
                            rna_file = NULL,
                            sample_name = NULL,
                            gene_list = NULL){
  mapping_type <- match.arg(mapping_type,c("genome","transcriptome"))

  # ============================================================================
  # extract data
  # ============================================================================
  if(mapping_type == "genome"){
    plyr::ldply(1:length(ribo_file),function(x){
      # ribo denisty
      # ribo_tmp <- data.table::fread(paste('2.density-data/',ribo_file[x],sep = ''),sep = "\t")
      ribo_tmp <- vroom::vroom(paste('2.density-data/',ribo_file[x],sep = ''),
                               delim = "\t",show_col_types = FALSE,
                               col_names = c('gene_name',"trans_id",'transpos','density'))

      # add colnames
      # colnames(ribo_tmp) <- c('gene_name',"trans_id",'transpos','density')
      # filter
      ribo_tmp <- ribo_tmp %>% dplyr::filter(gene_name %in% gene_list)
      # add type
      ribo_tmp$type <- 'ribo'

      if(!is.null(rna_file)){
        # rna coverage
        # rna_tmp <- data.table::fread(paste('2.density-data/',rna_file[x],sep = ''),sep = "\t")

        rna_tmp <- vroom::vroom(paste('2.density-data/',rna_file[x],sep = ''),
                                delim = "\t",show_col_types = FALSE,
                                col_names = c('gene_name',"trans_id",'transpos','density'))
        # add colnames
        # colnames(rna_tmp) <- c('gene_name',"trans_id",'transpos','density')
        # filter
        rna_tmp <- rna_tmp %>% dplyr::filter(gene_name %in% gene_list)
        # add type
        rna_tmp$type <- 'rna'

        # merge
        mer <- rbind(ribo_tmp,rna_tmp)
      }else{
        mer <- ribo_tmp
      }
      mer$sample <- sample_name[x]
      return(mer)
    }) -> df_gene
  }else if(mapping_type == "transcriptome"){
    plyr::ldply(1:length(ribo_file),function(x){
      # ribo denisty
      # ribo_tmp <- data.table::fread(paste('2.density-data/',ribo_file[x],sep = ''),
      #                               sep = "\t")[,c(1,2,4)]
      ribo_tmp <- vroom::vroom(paste('2.density-data/',ribo_file[x],sep = ''),
                               delim = "\t",show_col_types = FALSE,
                               col_names = c('id',"transpos",'none','density'),
                               col_select = c('id','transpos','density')) |>
        data.table::as.data.table()

      # add colnames
      # colnames(ribo_tmp) <- c('id','transpos','density')
      ribo_tmp[,c("gene_name") := data.table::tstrsplit(id,"|",fixed = TRUE)[1]]
      ribo_tmp[,c("trans_id") := data.table::tstrsplit(id,"|",fixed = TRUE)[3]]
      ribo_tmp <- ribo_tmp[,c("gene_name","trans_id",'transpos','density')]

      # filter
      ribo_tmp <- ribo_tmp %>% dplyr::filter(gene_name %in% gene_list)
      # add type
      ribo_tmp$type <- 'ribo'

      if(!is.null(rna_file)){
        # rna coverage
        # rna_tmp <- data.table::fread(paste('2.density-data/',rna_file[x],sep = ''),
        #                              sep = "\t")[,c(1,2,4)]
        rna_tmp <- vroom::vroom(paste('2.density-data/',rna_file[x],sep = ''),
                                delim = "\t",show_col_types = FALSE,
                                col_names = c('id',"transpos",'none','density'),
                                col_select = c('id','transpos','density')) |>
          data.table::as.data.table()

        # add colnames
        # colnames(rna_tmp) <- c('id','transpos','density')
        rna_tmp[,c("gene_name") := data.table::tstrsplit(id,"|",fixed = TRUE)[1]]
        rna_tmp[,c("trans_id") := data.table::tstrsplit(id,"|",fixed = TRUE)[3]]
        rna_tmp <- rna_tmp[,c("gene_name","trans_id",'transpos','density')]

        # filter
        rna_tmp <- rna_tmp %>% dplyr::filter(gene_name %in% gene_list)
        # add type
        rna_tmp$type <- 'rna'

        # merge
        mer <- rbind(ribo_tmp,rna_tmp)
      }else{
        mer <- ribo_tmp
      }
      mer$sample <- sample_name[x]
      return(mer)
    }) -> df_gene
  }
}




#' Preprocess SAM files to calculate count and TPM data
#'
#' This function takes SAM files as input, processes it using JuliaCall package
#' and generates count and TPM data.
#'
#' @param sam_file A character vector of SAM file paths.
#' @param out_file A character vector of output file names.
#' @param type Specifies the type of the input file. Default is 'ribo'.
#' @return NULL
#' @examples
#' \dontrun{
#' pre_count_tpm_data(sam_file = c("path/to/samfile1.sam","path/to/samfile2.sam"),
#'                    out_file = c("output1.txt","output2.txt"),
#'                    type = "rna")
#' }
#'
#' @export
pre_count_tpm_data <- function(sam_file = NULL,
                               out_file = NULL,
                               type = c("ribo","rna")){
  type <- match.arg(type,c("ribo","rna"))

  if(!dir.exists("4.expression-data")){
    dir.create("4.expression-data")
  }

  JuliaCall::julia_setup(installJulia = TRUE)

  JuliaCall::julia_library("XAM")

  script_path <- paste0('include("',
                        system.file("extdata", "CalculateCountTPM.jl",
                                    package = "RiboProfiler"),
                        '")',collapse = "")

  # choose function
  CalculateCountTPM <- JuliaCall::julia_eval(script_path)

  # excute function
  lapply(seq_along(sam_file), function(x){
    outFile_tmp = paste("4.expression-data/",out_file[x],sep = "")
    CalculateCountTPM(inputFile = sam_file[x],
                      outputFile = outFile_tmp,
                      inputType = type)
    message(paste(sam_file[x]," has been processed!",sep = ""))
  }) -> tmp

  return(NULL)
}


#' Load expression data from text files
#'
#' This function loads gene expression data from text files in a specified
#' directory. It returns both count and TPM (transcripts per million) matrices
#' for all samples.
#'
#' @param sample_name A character vector containing sample names. If NULL, the
#' function will use the file names as sample names.
#'
#' @return A list containing two matrices:
#' \describe{
#'   \item{count_matrix}{A matrix of gene expression counts for all samples.}
#'   \item{tpm_matrix}{A matrix of gene expression TPM values for all samples.}
#' }
#'
#' @export
load_expression_data <- function(sample_name = NULL){
  # load data
  file <- list.files('4.expression-data/','.txt')
  message("Expression input files: ")
  message(paste0(file,sep = "\n"))

  if(is.null(sample_name)){
    sample_name <- sapply(strsplit(file,split = '\\.'),'[',1)
  }else{
    sample_name <- sample_name
  }

  # ============================================================================
  # extract count data
  lapply(1:length(file),function(x){
    tmp <- data.table::fread(paste('4.expression-data/',file[x],sep = ''))[,-3]
    # add sample
    sample <- sample_name[x]
    colnames(tmp) <- c("gene_name",sample)

    return(tmp)
  }) -> count_list

  # merge
  all_count <- Reduce(function(x,y,...){merge(x,y,by = "gene_name",all = TRUE,...)},count_list)
  all_count[is.na(all_count)] <- 0

  # ============================================================================
  # extract tpm data
  lapply(1:length(file),function(x){
    tmp <- data.table::fread(paste('4.expression-data/',file[x],sep = ''))[,-2]
    # add sample
    sample <- sample_name[x]
    colnames(tmp) <- c("gene_name",sample)

    return(tmp)
  }) -> tpm_list

  # merge
  all_tpm <- Reduce(function(x,y,...){merge(x,y,by = "gene_name",all = TRUE,...)},tpm_list)
  all_tpm[is.na(all_tpm)] <- 0

  # return
  return(list(count_matrix = all_count,
              tpm_matrix = all_tpm))
}
