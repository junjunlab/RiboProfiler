#' Fetch tRNA and rRNA sequences from NCBI database
#'
#' This function allows users to fetch tRNA and rRNA sequences for a given
#' species from the NCBI nucleotide database.
#'
#' @param species A character string specifying the scientific name of the target
#' species. Defaults to "Homo sapiens".
#' @param output_file A character string specifying the name of the output file.
#' If not provided, the default name is "<species>_trRNA.fa".
#'
#' @return A FASTA-formatted file containing the tRNA and rRNA sequences for the
#' specified species is created.
#'
#' @references
#' \url{https://www.ncbi.nlm.nih.gov/}
#'
#' @examples
#' \dontrun{
#' # Fetch tRNA and rRNA sequences for human
#' fetch_trRNA_from_NCBI("Homo sapiens")
#'
#' # Fetch tRNA and rRNA sequences for mouse and save to a custom file name
#' fetch_trRNA_from_NCBI("Mus musculus", "my_mouse_trRNA.fa")
#' }
#'
#' @export
fetch_trRNA_from_NCBI <- function(species = "Homo sapiens",
                                  output_file = NULL){
  # search
  rs <- rentrez::entrez_search(db = "nucleotide",
                               term = paste("(",species,"[ORGN] AND rRNA[FILT]) OR (",
                                            species,"[ORGN] AND tRNA[FILT])",sep = ""),
                               retmax = 20000,
                               use_history = TRUE)

  # download sequence
  if(is.null(output_file)){
    output_file <- paste(gsub(" ",replacement = "_",species),"_trRNA.fa",sep = "")
    file.create(output_file)
  }else{
    file.create(output_file)
  }

  # x = 2
  ids = rs$ids
  fa <- rentrez::entrez_fetch(db = "nucleotide",rettype = "fasta",
                              web_history = rs$web_history)
  write(fa,file = output_file,append = TRUE)

  raw_fa <- Biostrings::readDNAStringSet(output_file)
  Biostrings::writeXStringSet(raw_fa,filepath = output_file)

  cli::cat_boxx("Download tRNA and rRNA sequences from NCBI has finished!",
                col = "red")

  cli::cat_bullet("Here are some common species Latin names:\n
    人类(Homo sapiens), 大鼠(Rattus norvegicus), 小鼠(Mus musculus),
    斑马鱼(Danio rerio), 红毛猩猩(Pan troglodytes), 家犬(Canis familiaris),
    草履虫(Dictyostelium discoideum), 猴子(Macaca mulatta), 红松鼠(Tamiasciurus hudsonicus),
    家兔(Oryctolagus cuniculus), 黄鼠狼(Cricetulus griseus), 南方豹猫(Prionailurus rubiginosus),
    畜牛(Bos taurus), 绿猴(Chlorocebus sabaeus), 绵羊(Ovis aries),
    猪(Sus scrofa), 验血鸟(Taeniopygia guttata), 萨摩耶犬(Canis lupus familiaris),
    膜骨鱼(Takifugu rubripes), 猫(Felis catus), 银狐(Vulpes vulpes),
    水稻(Oryza sativa), 吸血蝙蝠(Desmodus rotundus), 巴西三带蚊(Aedes aegypti),
    印度大象(Elephas maximus), 狼(Canis lupus), 仓鼠(Cricetulus griseus),
    刺参(Strongylocentrotus purpuratus), 山羊(Capra hircus), 兔子(Oryctolagus cuniculus),
    黄猴(Macaca fascicularis), 石斑鱼(Epinephelus coioides), 柴犬(Canis lupus familiaris),
    蚯蚓(Eisenia fetida), 小萤火虫(Luciola cruciata), 白化病毒(White spot syndrome virus),
    河北地蜂(Apis cerana), 喜马拉雅兔(Ochotona curzoniae), 裸鼠(Heterocephalus glaber),
    马(Equus caballus), 胡萝卜野生种(Daucus carota subsp. carota), 水牛(Bubalus bubalis),
    歌鸲(Erithacus rubecula), 美国黑熊(Ursus americanus), 鲨鱼(Callorhinchus milii),
    大黄蜂(Vespa mandarinia), 蜡螟(Galleria mellonella), 黄斑蝶(Papilio xuthus),
    皇家企鹅(Aptenodytes forsteri), 银河野猪(Sus scrofa)",bullet = "play",bullet_col = "purple",
    background_col = "grey98",col = "#FC4F00")
}



#' Build an index for tRNA and rRNA sequences using Bowtie2.
#'
#' This function builds an index for tRNA and rRNA sequences using Bowtie2.
#' The resulting index can be used for mapping reads against tRNA/rRNA sequences
#' to remove them from downstream analysis, as they can interfere with gene
#' expression measurement.
#'
#' @param trRNA_file A character string specifying the path to the FASTA file
#' containing tRNA and rRNA sequences.
#' @param prefix A character string specifying the prefix for the output index
#' files.
#' @param threads An integer specifying the number of threads to use for building
#' the index.
#'
#' @return Nothing is returned; the function saves the index files in the
#' "0.index-data/rtRNA-index/" directory.
#'
#' @examples
#' \dontrun{
#' trRNA_index_build("path/to/trRNA.fasta", "trRNA_bt2_index", 8)
#' }
#'
#' @export
trRNA_index_build <- function(trRNA_file = NULL,prefix = NULL,
                              threads = 8){
  # create directory
  if(!dir.exists("0.index-data/rtRNA-index/")){
    dir.create("0.index-data/rtRNA-index/",recursive = TRUE)
  }

  # Rbowtie2::bowtie2_build(references = trRNA_file,
  #                         bt2Index = paste("0.index-data/rtRNA-index/",prefix,sep = ""),
  #                         "--threads",threads,"--quiet",overwrite = TRUE)

  # download bowtie2 for windows
  bowtie2 = "https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-mingw-x86_64.zip"

  if(!dir.exists("bowtie2-2.5.1-mingw-x86_64")){
    if(suppressWarnings(download.file(url = bowtie2,destfile = "./"))){
      cli::cli_alert_warning(paste("Please download bowtie2 mannually with links:",
                                   bowtie2,
                                   " and unzip file in current directory.",
                                   sep = " "))
    }else{
      unzip("./bowtie2-2.5.1-mingw-x86_64.zip")
    }
  }

  bowtie_path = "python ./bowtie2-2.5.1-mingw-x86_64/bowtie2-build"

  # build index
  # system2(bowtie_path,c(paste("--threads ",threads,sep = ""),
  #                       "--quiet",trRNA_file,
  #                       paste("0.index-data/rtRNA-index/",prefix,sep = "")))

  system(paste(bowtie_path,
               paste("--threads ",threads,sep = ""),
               "--quiet",trRNA_file,
               paste("0.index-data/rtRNA-index/",prefix,sep = ""),
               sep = " "))

  cli::cat_boxx("Building index for tRNA and rRNA has finished!",
                col = "red")
}



#' Build HISAT2 index for reference genome
#'
#' This function builds a HISAT2 index for the given reference genome using the
#' HISAT2 tool.
#'
#' @param reference_file path to the reference genome file in FASTA format
#' @param prefix prefix for the output files
#' @param threads number of CPU threads to use for building the index
#' @param hisat2_build_params a list of additional parameters to pass to HISAT2
#' build command
#'
#' @examples
#' \dontrun{
#' # Build a HISAT2 index for human reference genome
#' reference_index_build(reference_file = "hg38.fa", prefix = "hg38")
#' }
#'
#' @export
reference_index_build <- function(reference_file = NULL,prefix = NULL,
                                  threads = 8,hisat2_build_params = list()){
  # create directory
  if(!dir.exists("0.index-data/ref-index/")){
    dir.create("0.index-data/ref-index/",recursive = TRUE)
  }

  # build index
  do.call(Rhisat2::hisat2_build,
          modifyList(list(references = reference_file,
                          outdir = "0.index-data/ref-index/",
                          prefix = prefix,
                          p = threads,quiet = TRUE,force = TRUE),
                     hisat2_build_params))

  cli::cat_boxx("Building index for reference has finished!",
                col = "red")
}


#' Batch Adapter Removal using Rfastp
#'
#' This function removes adapters from a batch of paired-end fastq files using Rfastp.
#'
#' @param fastq_file1 A vector of file paths to the first read in each pair.
#' @param fastq_file2 A vector of file paths to the second read in each pair.
#' Default is NULL, indicating single-end reads.
#' @param output_dir The directory to save the output files in.
#' @param output_name A vector of names for the output files. Must have the same
#' length as fastq_file1 and fastq_file2.
#' @param fastp_params A list of parameters to pass to Rfastp.
#' See https://github.com/lh3/fastp for more information.
#'
#' @return A list of QC summaries for each input file.
#'
#' @examples
#' \dontrun{
#' # Batch adapter removal for paired-end reads
#' batch_adapter_remove(fastq_file1 = c("read1_1.fastq", "read1_2.fastq"),
#'                      fastq_file2 = c("read2_1.fastq", "read2_2.fastq"),
#'                      output_dir = "output/",
#'                      output_name = c("out1.fastq", "out2.fastq"))
#'
#' # Batch adapter removal for single-end reads
#' batch_adapter_remove(fastq_file1 = c("read1.fastq", "read2.fastq"),
#'                      output_dir = "output/",
#'                      output_name = c("out1.fastq", "out2.fastq"))
#' }
#'
#' @export
batch_adpator_remove <- function(fastq_file1 = NULL,
                                 fastq_file2 = NULL,
                                 output_dir = NULL,
                                 output_name = NULL,
                                 fastp_params = list()){
  # check read 2
  if(is.null(fastq_file2)){
    read2 <- rep("",length(fastq_file1))
  }else{
    read2 <- fastq_file2
  }

  # x = 1
  lapply(seq_along(fastq_file1), function(x){
    json_report <- do.call(Rfastp::rfastp,
                           modifyList(list(read1 = fastq_file1[x],
                                           read2 = read2[x],
                                           outputFastq = paste(output_dir,output_name[x],sep = ""),
                                           verbose = FALSE),
                                      fastp_params))

    cli::cat_bullet(paste(output_name[x],"has been processed!",sep = " "),
                    bullet = "play",bullet_col = "purple",
                    background_col = "grey98",col = "#FC4F00")
    return(json_report)
  }) -> qc_summary_list
  names(qc_summary_list) <- output_name
}



#' Bowtie2 Alignment
#'
#' This function performs Bowtie2 alignment of sequencing reads to a reference
#' genome. It requires the path to the Bowtie2 index files, input FASTQ files,
#' output file name, number of threads, and Bowtie2 parameters as inputs. If the
#' input FASTQ files are compressed, it automatically decompresses them before
#' alignment. The output includes both the mapping result in SAM format and the
#' mapping information in a separate text file.
#'
#' @param index Path to the Bowtie2 index files for the reference genome.
#' @param fq_file1 Path to the input FASTQ file containing the first set of
#' paired-end reads or single-end reads.
#' @param fq_file2 Path to the input FASTQ file containing the second set of
#' paired-end reads. Set to NULL if analyzing single-end reads.
#' @param output_file Path for the output SAM file. An additional text file
#' containing mapping information will also be generated with the same file name
#' but with "_mapinfo" appended to the end.
#' @param threads Number of CPU threads to use for alignment.
#' @param bowtie2_params Additional Bowtie2 parameters to be included in the command line.
#'
#' @return Prints a message indicating that the output file has been generated.
#' @export
bowtie2_align <- function(index = NULL,fq_file1 = NULL,fq_file2 = NULL,
                          output_file = NULL,threads = 1,bowtie2_params = NULL){
  # ============================================================================
  # download bowtie2
  bowtie2 = "https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-mingw-x86_64.zip"

  if(!dir.exists("bowtie2-2.5.1-mingw-x86_64")){
    if(suppressWarnings(download.file(url = bowtie2,destfile = "./"))){
      cli::cli_alert_warning(paste("Please download bowtie2 mannually with links:",
                                   bowtie2,
                                   " and unzip file in current directory.",
                                   sep = " "))
    }else{
      unzip("./bowtie2-2.5.1-mingw-x86_64.zip")
    }
  }

  bowtie_path = "./bowtie2-2.5.1-mingw-x86_64/bowtie2"

  # ============================================================================
  tmp_params = c(paste("-x ",index,sep = ""),
                 paste("--threads ",threads,sep = ""),
                 paste("-S ",output_file,sep = ""),
                 bowtie2_params)

  # decompress fsatq file
  if(is.null(fq_file2)){
    tmp_file = sapply(strsplit(fq_file1,split = ".gz|.gzip"),"[",1)

    if(!file.exists(tmp_file)){
      R.utils::gunzip(fq_file1,remove = FALSE)
    }

    final_params = c(paste("-U ",tmp_file,sep = ""),tmp_params)
  }else{
    tmp_file1 = sapply(strsplit(fq_file1,split = ".gz|.gzip"),"[",1)
    tmp_file2 = sapply(strsplit(fq_file2,split = ".gz|.gzip"),"[",1)

    if(!file.exists(tmp_file1)){
      R.utils::gunzip(fq_file1,remove = FALSE)
    }

    if(!file.exists(tmp_file2)){
      R.utils::gunzip(fq_file2,remove = FALSE)
    }

    final_params = c(paste("-1 ",tmp_file1,sep = ""),
                     paste("-2 ",tmp_file2,sep = ""),
                     tmp_params)
  }

  # mapping
  map <- system2(bowtie_path,final_params,
                 stdout = TRUE)

  # output map info
  map_info <- paste0(map,collapse = "\n")
  write(map_info,file = paste(output_file,"_mapinfo.txt",sep = ""))

  # remove fastq
  if(is.null(fq_file2)){
    file.remove(tmp_file)
  }else{
    file.remove(tmp_file1)
    file.remove(tmp_file2)
  }

  # print
  cli::cat_bullet(paste(output_file,"has been processed!",sep = " "),
                  bullet = "play",bullet_col = "purple",
                  background_col = "grey98",col = "#FC4F00")
}



#' Batch align FASTQ files using HISAT2
#'
#' This function performs batch alignment of FASTQ files using HISAT2.
#'
#' @param index A character string specifying the path to the HISAT2 index.
#' @param fq_file1 A character vector specifying the paths to the first input
#' FASTQ files.
#' @param fq_file2 A character vector specifying the paths to the second input
#' FASTQ files (optional).
#' @param output_dir A character string specifying the path to the output directory.
#' @param output_file A character vector specifying the names for the output files.
#' @param threads An integer specifying the number of threads to use in the
#' alignment process (default: 1).
#' @param hisat2_params A list specifying additional parameters for HISAT2
#' (default: empty list).
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Align a single-end FASTQ file
#' batch_hisat_align(index = "path/to/index",
#'                   fq_file1 = "path/to/fq_file.fastq.gz",
#'                   output_dir = "path/to/output/",
#'                   output_file = "aligned.sam")
#'
#' # Align paired-end FASTQ files
#' batch_hisat_align(index = "path/to/index",
#'                   fq_file1 = c("path/to/fq_file1_1.fastq.gz", "path/to/fq_file2_1.fastq.gz"),
#'                   fq_file2 = c("path/to/fq_file1_2.fastq.gz", "path/to/fq_file2_2.fastq.gz"),
#'                   output_dir = "path/to/output/",
#'                   output_file = c("aligned_1.sam", "aligned_2.sam"))
#' }
#'
#' @export
batch_hisat_align <- function(index = NULL,
                              fq_file1 = NULL,fq_file2 = NULL,
                              output_dir = NULL,
                              output_file = NULL,threads = 1,hisat2_params = list()){
  lapply(seq_along(fq_file1), function(x){
    # check whether supply paired reads
    if(is.null(fq_file2)){
      if(endsWith(fq_file1[x],".gz|.gzip")){
        R.utils::gunzip(fq_file1[x],remove = FALSE)
        sequences <- sapply(strsplit(fq_file1[x],split = ".gz|.gzip"),"[",1)
      }else{
        sequences <- fq_file1[x]
      }
    }else{
      # ================================================
      # paired end
      if(endsWith(fq_file1[x],".gz|.gzip")){
        R.utils::gunzip(fq_file1[x],remove = FALSE)
        sequences_1 <- sapply(strsplit(fq_file1[x],split = ".gz|.gzip"),"[",1)
      }else{
        sequences_1 <- fq_file1[x]
      }

      if(endsWith(fq_file2[x],".gz|.gzip")){
        R.utils::gunzip(fq_file2[x],remove = FALSE)
        sequences_2 <- sapply(strsplit(fq_file2[x],split = ".gz|.gzip"),"[",1)
      }else{
        sequences_2 <- fq_file2[x]
      }

      sequences <- list(sequences_1,sequences_2)
    }

    outfile <- paste(output_dir,output_file[x],sep = "")

    # mapping now
    tmp <- do.call(Rhisat2::hisat2,
                   modifyList(list(sequences = sequences,
                                   index = index,
                                   p = threads,force = TRUE,
                                   outfile = paste(outfile,".sam",sep = ""),
                                   "summary-file" = paste(outfile,"_mapinfo.txt",sep = "")),
                              hisat2_params))

    # print
    cli::cat_bullet(paste(output_file[x],"has been processed!",sep = " "),
                    bullet = "play",bullet_col = "orange",
                    background_col = "grey98",col = "#00235B")
  }) -> tmp
}



#' Convert BAM files to bigWig format
#'
#' This function takes a set of BAM files and converts them to bigWig format.
#' The BAM files are read in using the \code{GenomicAlignments} package, extended
#' to an appropriate fragment length, converted to a GRanges object, and coverage
#' is calculated. The coverage data is then normalized to reads per million (RPM)
#' and exported as a bigWig file using the \code{rtracklayer} package.
#'
#' @param bam_file A character vector giving the paths to the input BAM files.
#' @param bw_file A character vector giving the paths to the output bigWig files.
#' @param paired The sequencing type of reads, defaults False(single-end).
#'
#' @details The function uses the "readGAlignments" function from the
#' \code{GenomicAlignments} package to read in each BAM file. The resulting
#' object is then converted to a GRanges object using the "as" function. Coverage
#' is obtained using the "coverage" function from the \code{IRanges} package.
#' The coverage data is then normalized to RPM by dividing by the total number of
#' reads and multiplying by one million. The normalized data is stored as a
#' SimpleRleList object, which is then exported to a bigWig file using the
#' \code{export.bw} function from the \code{rtracklayer} package.
#'
#' @return The function invisibly returns a list containing the output file paths.
#'
#' @import GenomicAlignments
#' @importFrom IRanges coverage
#' @importFrom rtracklayer export.bw
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#' # Convert a single BAM file to bigWig format
#' batch_bam2bigwig(bam_file = "/path/to/my.bam",
#'                  bw_file = "/path/to/my.bw")
#'
#' # Convert multiple BAM files to bigWig format
#' batch_bam2bigwig(bam_file = c("/path/to/my1.bam", "/path/to/my2.bam"),
#'                  bw_file = c("/path/to/my1.bw", "/path/to/my2.bw"))
#' }
#'
#' @export
batch_bam2bigwig <- function(bam_file = NULL,
                             bw_file = NULL,
                             paired = FALSE){
  # ============================================================================
  # define transform function
  # ============================================================================
  # create a function to read in data
  # extend reads, get coverage, normalize to RPM
  # and export to bw
  bam2bw <- function(bamfile,bwfile,paired){

    # cat("opening:", bamfile, sep="\n")
    if (paired) {
      bd <- GenomicAlignments::readGAlignmentPairs(bamfile)
    } else {
      bd <- GenomicAlignments::readGAlignments(bamfile)
    }


    # cat("convert to GRanges\n")
    mygr <- methods::as(bd,"GRanges")
    totalReads <- length(mygr)

    # cat("getting coverage\n")
    # get coverage
    cov <- IRanges::coverage(mygr)

    # rpm
    rpm <- lapply(cov, function(x) signif(10^6 * x/totalReads,3))
    rpm <- as(rpm,"SimpleRleList")

    # export rpm to bigWig
    rtracklayer::export.bw(rpm, bwfile)
  }


  # for each element of our vector, call the bam2bw function
  lapply(seq_along(bam_file), function(x){
    bam2bw(bam_file[x],bw_file[x],paired)

    # print
    cli::cat_bullet(paste(bam_file[x],"has been processed!",sep = " "),
                    bullet = "play",bullet_col = "orange",
                    background_col = "grey98",col = "#00235B")
  }) -> tmp
}



#' Convert SAM files to BAM format
#'
#' This function converts one or more SAM files to BAM format using
#' Rsamtools::asBam function.
#'
#' @param sam_file Character vector containing the path(s) to the input SAM file(s).
#' @param bam_file Character vector containing the path(s) to where the output
#' BAM file(s) will be saved.
#' @return None
#'
#' @examples
#' \dontrun{
#' # Convert a single SAM file to BAM format
#' batch_sam2bam("my_sam_file.sam", "my_bam_file")
#'
#' # Convert multiple SAM files to BAM format
#' batch_sam2bam(c("sam_file1.sam", "sam_file2.sam"), c("bam_file1", "bam_file2"))
#' }
#'
#' @export
batch_sam2bam <- function(sam_file = NULL,
                          bam_file = NULL){
  lapply(seq_along(sam_file), function(x){
    Rsamtools::asBam(file = sam_file[x],destination = bam_file[x],overwrite = TRUE)

    # print
    cli::cat_bullet(paste(sam_file[x],"has been processed!",sep = " "),
                    bullet = "play",bullet_col = "orange",
                    background_col = "grey98",col = "#00235B")
  }) -> tmp
}



#' Calculate TPM from feature counts object
#'
#' This function takes a feature counts object and calculates the transcripts per
#' million (TPM) for each gene.
#'
#' @param fc_obj Feature counts object containing count data and gene annotation
#' which produced from subread featureCounts function.
#'
#' @return A data frame containing gene IDs, gene names, gene biotypes, and TPM values.
#'
#' @examples
#'
#' \dontrun{
#' # Calculate TPM using count_to_tpm()
#' tpm_df <- count_to_tpm(example_fc_obj)
#' }
#'
#' @export
count_to_tpm <- function(fc_obj = NULL){
  # counts
  my_counts <- fc_obj$counts

  # annotation info
  gene_anno <- fc_obj$annotation[,-2:-6]
  colnames(gene_anno) <- c("gene_id","gene_name","gene_biotype")

  # calculate tpm
  kb <- fc_obj$annotation$Length/1000
  rpk <- my_counts/kb
  tpm <- t(t(rpk)/colSums(rpk)*1000000)
  tpm <- as.data.frame(tpm)
  tpm$gene_id <- rownames(tpm)
  tpm_anno <- merge(tpm,gene_anno,by = 'gene_id')

  return(tpm_anno)
}
