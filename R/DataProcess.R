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


#' pre_qc_data function
#'
#' This function performs quality control (QC) analysis on ribosome profiling
#' data. The input is a SAM/BAM file generated from ribosome profiling data and the
#' output is a QC result in a text file format.
#'
#' @param mapping_type The mapping type for your sam files, "genome" or "transcriptome".
#' @param julia_path The julia program path on the computer.
#' @param longest_trans_file A string specifying the path to the longest transcript
#' file.
#' @param sam_file A character vector specifying the paths to the SAM files.
#' @param bam_file A character vector specifying the paths to the BAM files.
#' @param out_file A character vector specifying the paths to the output QC result
#' files.
#' @param seq_type The sequencing type for fastq files, "singleEnd" or "pairedEnd".
#'
#' @return The function does not return a value, but outputs QC results in text
#' files.
#'
#' @importFrom JuliaCall julia_eval
#'
#' @examples
#' \dontrun{
#' pre_qc_data(longest_trans_file = "path/to/longestTransInfo.txt",
#'             sam_file = c("path/to/sample1.sam", "path/to/sample2.sam"),
#'             out_file = c("path/to/sample1_QC_result.txt", "path/to/sample2_QC_result.txt"))
#' }
#'
#' @export
pre_qc_data <- function(mapping_type = c("genome","transcriptome"),
                        julia_path = NULL,
                        longest_trans_file = NULL,
                        sam_file = NULL,
                        bam_file = NULL,
                        out_file = NULL,
                        seq_type = c("pairedEnd","singleEnd")){
  mapping_type <- match.arg(mapping_type,c("genome","transcriptome"))
  seq_type <- match.arg(seq_type,c("pairedEnd","singleEnd"))

  if(!dir.exists("1.QC-data")){
    dir.create("1.QC-data")
  }

  # check julia
  if(!is.null(julia_path)){
    options(JULIA_HOME = julia_path)
  }

  options(timeout = max(1000, getOption("timeout")))
  JuliaCall::julia_setup(installJulia = TRUE)
  JuliaCall::julia_install_package_if_needed("XAM")
  JuliaCall::julia_library("XAM")

  # check input file type
  if(!is.null(sam_file) & is.null(bam_file)){
    script_path <- paste0('include("',
                          system.file("extdata", "prepareQCdata.jl",
                                      package = "RiboProfiler"),
                          '")',collapse = "")
  }else{
    script_path <- paste0('include("',
                          system.file("extdata", "prepareQCdataForBam.jl",
                                      package = "RiboProfiler"),
                          '")',collapse = "")
  }


  # choose function
  JuliaCall::julia_eval(script_path)

  if(mapping_type == "genome"){
    prepareQCdata <- JuliaCall::julia_eval("prepareQCdata")

    # inputfile
    if(!is.null(sam_file) & is.null(bam_file)){
      inFile = paste0(sam_file,collapse = ",")
    }else{
      inFile = paste0(bam_file,collapse = ",")
    }

    # excute function
    outFile_tmp = paste("1.QC-data/",out_file,sep = "")
    prepareQCdata(longestTransInfo = longest_trans_file,
                  inFile = inFile,
                  outFile = paste0(outFile_tmp,collapse = ","),
                  seqType = seq_type)
  }else if(mapping_type == "transcriptome"){
    prepareQCdata <- JuliaCall::julia_eval("prepareQCdata_ontrans")

    # inputfile
    if(!is.null(sam_file) & is.null(bam_file)){
      inFile = sam_file
    }else{
      inFile = bam_file
    }

    lapply(1:length(inFile), function(x){
      # excute function
      outFile_tmp = paste("1.QC-data/",out_file[x],sep = "")
      prepareQCdata(inFile = inFile[x],
                    outFile = outFile_tmp,
                    seqType = seq_type)
      message(paste(inFile[x]," has been processed!",sep = ""))
    }) -> tmp
  }

}


#' Calculate ribosome density data
#'
#' This function calculates the ribosome density data using XAM package in Julia.
#'
#' @param mapping_type The mapping type for your sam/bam files, "genome" or "transcriptome".
#' @param sam_file A character vector of SAM file paths.
#' @param bam_file A character vector of BAM file paths.
#' @param out_file A character vector of output file names.
#' @param min The minimum length of reads to be considered for calculating density
#' (default is 23).
#' @param max The maximum length of reads to be considered for calculating density
#' (default is 35).
#'
#' @return No explicit return value. Output files are written to "2.density-data" directory.
#'
#' @importFrom JuliaCall julia_setup julia_library julia_eval
#'
#' @examples
#' \dontrun{
#' # Calculate ribosome density data for a single SAM file
#' pre_ribo_density_data(sam_file = "sample.sam", out_file = "sample_density.txt")
#'
#' # Calculate ribosome density data for multiple SAM files
#' pre_ribo_density_data(sam_file = c("sample1.sam", "sample2.sam"),
#'                       out_file = c("sample1_density.txt", "sample2_density.txt"))
#' }
#'
#' @export
pre_ribo_density_data <- function(mapping_type = c("genome","transcriptome"),
                                  sam_file = NULL,
                                  bam_file = NULL,
                                  out_file = NULL,
                                  min = 23,max = 35){
  mapping_type <- match.arg(mapping_type,c("genome","transcriptome"))

  if(!dir.exists("2.density-data")){
    dir.create("2.density-data")
  }

  JuliaCall::julia_setup(installJulia = TRUE)

  JuliaCall::julia_library("XAM")

  # check input file type
  if(!is.null(sam_file) & is.null(bam_file)){
    script_path <- paste0('include("',
                          system.file("extdata", "CalculateRibosomeDensity.jl",
                                      package = "RiboProfiler"),
                          '")',collapse = "")
  }else{
    script_path <- paste0('include("',
                          system.file("extdata", "CalculateRibosomeDensityForBam.jl",
                                      package = "RiboProfiler"),
                          '")',collapse = "")
  }


  # choose function
  JuliaCall::julia_eval(script_path)

  if(mapping_type == "genome"){
    calculateRibosomeDensity <- JuliaCall::julia_eval("CalculateRibosomeDensity")
  }else if(mapping_type == "transcriptome"){
    calculateRibosomeDensity <- JuliaCall::julia_eval("CalculateRibosomeDensity_ontrans")
  }

  # excute function
  if(!is.null(sam_file) & is.null(bam_file)){
    inputFile <- sam_file
  }else{
    inputFile <- bam_file
  }

  # loop analysis
  lapply(seq_along(inputFile), function(x){
    outFile_tmp = paste("2.density-data/",out_file[x],sep = "")
    calculateRibosomeDensity(inputFile = inputFile[x],
                             outputFile = outFile_tmp,
                             min = min,
                             max = max)
    message(paste(inputFile[x]," has been processed!",sep = ""))
  }) -> tmp

  return(NULL)
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



#' load_qc_data function
#'
#' A function to load quality control (QC) data from files and create a merged
#' data frame.
#'
#' @param sample_name A vector of file names for the samples to be loaded.
#' @param group_name The name of the group to which the loaded samples belong.
#'
#' @return Returns a merged data frame containing QC data for all loaded samples.
#'
#' @importFrom plyr ldply
#' @importFrom data.table fread
#' @export
load_qc_data <- function(sample_name = NULL,
                         group_name = NULL){
  # load data
  file <- list.files('1.QC-data/','.qc.txt')

  message("QC input files: ")
  message(paste0(file,sep = "\n"))

  if(is.null(sample_name)){
    sample_name <- sapply(strsplit(file,split = '\\.'),'[',1)
  }else{
    sample_name <- sample_name
  }

  if(is.null(group_name)){
    group_name <- NA
  }else{
    group_name <- group_name
  }

  # loop read file
  plyr::ldply(1:length(file),function(x){
    # tmp <- data.table::fread(paste('1.QC-data/',file[x],sep = ''))
    tmp <- vroom::vroom(file = paste('1.QC-data/',file[x],sep = ''),col_names = F,show_col_types = FALSE)

    colnames(tmp) <- c('length','framest','relst','framesp','relsp',
                       'feature','trans_pos','trans_id','counts')
    # add sample
    tmp$sample <- sample_name[x]
    # add group
    tmp$group <- group_name[x]

    return(tmp)
  }) -> dfqc

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


#' Calculate Metagene Data
#'
#' This function calculates metagene data using ribosome density files and gene
#' annotation.
#'
#' @param mapping_type The mapping type for your sam files, "genome" or "transcriptome".
#' @param gene_anno A data.frame containing gene annotation information which is
#' needed when mapping_type = "genome".
#' @param density_file A character vector of file names for ribosome density data.
#' @param out_file A character vector of output file names for metagene data.
#' @param mode A character indicating the mode of analysis, "st" or "sp".
#' Default is "st".
#' @param type A character indicating the type of analysis, "nt" or "codon". Default is "codon".
#' @param cdslength An integer indicating the length of CDS. Default is 600.
#' @param expression An integer indicating the minimum expression value. Default is 30.
#' @param exclude An integer indicating the number of nucleotides to exclude from
#' both ends of the coding region. Default is 90.
#'
#' @return NULL
#'
#' @examples
#' # Assume that we have a gene annotation data frame called 'gene_annot' and two ribosome density
#' # files called 'file1.bam' and 'file2.bam' in '2.density-data/' directory. The following code
#' # will calculate the metagene data and save the results in '3.metagene-data/' directory with
#' # names 'metagene1.tsv' and 'metagene2.tsv':
#' #
#' # pre_metagene_data(gene_anno = gene_annot,
#' #                   density_file = c("file1.bam", "file2.bam"),
#' #                   out_file = c("metagene1.txt", "metagene2.txt"))
#'
#' @export
pre_metagene_data <- function(mapping_type = c("genome","transcriptome"),
                              gene_anno = NULL,
                              density_file = NULL,
                              out_file = NULL,
                              mode = c("st","sp"),
                              type = "codon",
                              cdslength = 600,
                              expression = 30,
                              exclude = 90){
  mapping_type <- match.arg(mapping_type,c("genome","transcriptome"))
  mode <- match.arg(mode,c("st","sp"))

  if(!dir.exists("3.metagene-data")){
    dir.create("3.metagene-data")
  }

  JuliaCall::julia_setup(installJulia = TRUE)

  JuliaCall::julia_install_package_if_needed("DataStructures")
  JuliaCall::julia_library("DataStructures")

  script_path <- paste0('include("',
                        system.file("extdata", "MetageneAnalysis.jl",
                                    package = "RiboProfiler"),
                        '")',collapse = "")

  # choose function
  JuliaCall::julia_eval(script_path)
  if(mapping_type == "genome"){
    MetageneAnalysis <- JuliaCall::julia_eval("MetageneAnalysis")

    # excute function
    lapply(seq_along(density_file), function(x){
      inputFile_tmp = paste("2.density-data/",density_file[x],sep = "")
      outFile_tmp = paste("3.metagene-data/",out_file[x],sep = "")
      MetageneAnalysis(geneInfo = gene_anno,
                       inputFile = inputFile_tmp,
                       outputFile = outFile_tmp,
                       mode = mode,
                       type = type,
                       cdslength = as.integer(cdslength),
                       expression = as.integer(expression),
                       exclude = as.integer(exclude))
      message(paste(inputFile_tmp," has been processed!"))
    }) -> tmp
  }else if(mapping_type == "transcriptome"){
    MetageneAnalysis <- JuliaCall::julia_eval("MetageneAnalysis_ontrans")

    # excute function
    lapply(seq_along(density_file), function(x){
      inputFile_tmp = paste("2.density-data/",density_file[x],sep = "")
      outFile_tmp = paste("3.metagene-data/",out_file[x],sep = "")
      MetageneAnalysis(inputFile = inputFile_tmp,
                       outputFile = outFile_tmp,
                       mode = mode,
                       type = type,
                       cdslength = as.integer(cdslength),
                       expression = as.integer(expression),
                       exclude = as.integer(exclude))
      message(paste(inputFile_tmp," has been processed!"))
    }) -> tmp
  }


  return(NULL)
}

#' Load Metagene Data
#'
#' This function loads metagene data from text files in the specified directory.
#'
#' @param sample_name A character vector specifying the names of the input files.
#' If NULL, all files in the directory will be used.
#' @param group_name A character vector specifying the group names for each set
#' of input files. If NULL, no grouping will be performed.
#'
#' @return A data frame containing the loaded metagene data with columns for
#' position, density, sample name, and group name.
#'
#' @examples
#' \dontrun{
#' load_metagene_data(sample_name = c("sample1.txt", "sample2.txt"))
#' }
#'
#' @export
load_metagene_data <- function(sample_name = NULL,
                               group_name = NULL){
  # load data
  file <- list.files('3.metagene-data/','.txt')
  message("MetaGene input files: ")
  message(paste0(file,sep = "\n"))

  if(is.null(sample_name)){
    sample_name <- sapply(strsplit(file,split = '\\.'),'[',1)
  }else{
    sample_name <- sample_name
  }

  if(is.null(group_name)){
    group_name <- NA
  }else{
    group_name <- group_name
  }

  plyr::ldply(1:length(file),function(x){
    tmp <- data.table::fread(paste('3.metagene-data/',file[x],sep = ''))
    colnames(tmp) <- c('pos','density')
    # add sample
    tmp$sample <- sample_name[x]
    # add group
    tmp$group <- group_name[x]
    return(tmp)
  }) -> dfmeta
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
