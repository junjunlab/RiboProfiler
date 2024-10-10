#' pre_qc_data function
#'
#' This function performs quality control (QC) analysis on ribosome profiling
#' data. The input is a SAM/BAM file generated from ribosome profiling data and the
#' output is a QC result in a text file format.
#'
#' @param mapping_type The mapping type for your sam files, "genome" or "transcriptome".
#' @param julia_path The julia program path on the computer.
#' @param XAM_version The XAM version you have installed, "3"(0.3.1) or "4"(0.4.0).
#' @param longest_trans_file A string specifying the path to the longest transcript
#' file.
#' @param sam_file A character vector specifying the paths to the SAM files.
#' @param bam_file A character vector specifying the paths to the BAM files.
#' @param out_file_dir A character vector specifying the paths to the output data.
#' @param seq_type The sequencing type for fastq files, "singleEnd" or "pairedEnd".
#' @param assignType The reads assignment type, 5'end or 3'end, choices c("end5","end3"). Default 5'end.
#' @param out_file_dir A character vector specifying the paths to the output data.
#' @param out_file_prefix A character string specifying the prefix for the output files generated during
#' the QC analysis.
#' @param rep_name The replicates name. Default NULL.
#' @param group_name The group name. Default NULL.
#' @param has_created_data whether you have created raw data and save time to renalysis.
#'  Defaluts FALSE.
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
                        XAM_version = c("3","4"),
                        longest_trans_file = NULL,
                        sam_file = NULL,
                        bam_file = NULL,
                        out_file_dir = NULL,
                        out_file_prefix = NULL,
                        rep_name = NULL,
                        group_name = NULL,
                        has_created_data = FALSE,
                        assignType = c("end5","end3"),
                        seq_type = c("pairedEnd","singleEnd")){
  mapping_type <- match.arg(mapping_type,c("genome","transcriptome"))
  seq_type <- match.arg(seq_type,c("pairedEnd","singleEnd"))
  assignType <- match.arg(assignType,c("end5","end3"))
  XAM_version <- match.arg(XAM_version,c("3","4"))

  if(has_created_data == FALSE){
    if(!dir.exists("1.raw-data")){
      dir.create("1.raw-data")
    }

    # create out_file_dir
    if(!is.null(out_file_dir)){
      dir <- paste("1.raw-data/",out_file_dir,sep = "")
      if(!dir.exists(dir)){
        dir.create(dir)
      }
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
      # check XAM version
      if(XAM_version == "3"){
        prepareQCdata <- JuliaCall::julia_eval("prepareQCdata")
      }else{
        prepareQCdata <- JuliaCall::julia_eval("prepareQCdata2")
      }

      # inputfile
      if(!is.null(sam_file) & is.null(bam_file)){
        inFile = paste0(sam_file,collapse = ",")
      }else{
        inFile = paste0(bam_file,collapse = ",")
      }

      # excute function
      outFile_tmp = paste("1.raw-data/",out_file_dir,out_file_prefix,".raw.txt",sep = "")
      prepareQCdata(longestTransInfo = longest_trans_file,
                    inFile = inFile,
                    assignType = assignType,
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
        outFile_tmp = paste("1.raw-data/",out_file_dir,out_file_prefix[x],".raw.txt",sep = "")
        prepareQCdata(inFile = inFile[x],
                      outFile = outFile_tmp,
                      assignType = assignType,
                      seqType = seq_type)
        message(paste(inFile[x]," has been processed!",sep = ""))
      }) -> tmp
    }
  }


  # ============================================================================
  # load raw data in R
  # ============================================================================
  # loop read file
  file <- paste("1.raw-data/",out_file_dir,out_file_prefix,".raw.txt",sep = "")
  plyr::ldply(seq_along(file),function(x){
    tmp <- vroom::vroom(file = file[x],col_names = F,show_col_types = FALSE)

    colnames(tmp) <- c('len','framest','relst','framesp','relsp',
                       'feature','trans_pos','trans_id','counts')

    # add sample
    tmp$sample <- out_file_prefix[x]

    if(is.null(rep_name)){
      tmp$rep <- out_file_prefix[x]
    }else{
      tmp$rep <- rep_name[x]
    }

    if(is.null(group_name)){
      tmp$group <- NA
    }else{
      tmp$group <- group_name[x]
    }

    return(tmp)
  }) -> dfqc

}






#' #' load_qc_data function
#' #'
#' #' A function to load quality control (QC) data from files and create a merged
#' #' data frame.
#' #'
#' #' @param sample_name A vector of file names for the samples to be loaded.
#' #' @param group_name The name of the group to which the loaded samples belong.
#' #'
#' #' @return Returns a merged data frame containing QC data for all loaded samples.
#' #'
#' #' @importFrom plyr ldply
#' #' @importFrom data.table fread
#' #' @export
#' load_qc_data <- function(sample_name = NULL,
#'                          group_name = NULL){
#'   # load data
#'   file <- list.files('1.raw-data/','.raw.txt')
#'
#'   message("QC input files: ")
#'   message(paste0(file,sep = "\n"))
#'
#'   if(is.null(sample_name)){
#'     sample_name <- sapply(strsplit(file,split = '\\.raw.txt'),'[',1)
#'   }else{
#'     sample_name <- sample_name
#'   }
#'
#'   if(is.null(group_name)){
#'     group_name <- NA
#'   }else{
#'     group_name <- group_name
#'   }
#'
#'   # loop read file
#'   plyr::ldply(1:length(file),function(x){
#'     # tmp <- data.table::fread(paste('1.raw-data/',file[x],sep = ''))
#'     tmp <- vroom::vroom(file = paste('1.raw-data/',file[x],sep = ''),col_names = F,show_col_types = FALSE)
#'
#'     colnames(tmp) <- c('length','framest','relst','framesp','relsp',
#'                        'feature','trans_pos','trans_id','counts')
#'     # add sample
#'     tmp$sample <- sample_name[x]
#'     # add group
#'     tmp$group <- group_name[x]
#'
#'     return(tmp)
#'   }) -> dfqc
#'
#' }
