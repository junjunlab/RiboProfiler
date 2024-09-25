#' FetchSequence
#'
#' @param object ribosomeObj object.
#' @param output_prefix A character string specifying the path to the output file.
#' @param type A character string specifying the type of sequence to fetch ('cds' or 'exon').
#' @param coding_type A character string specifying the coding type ('NT' or 'AA').
#' @param pythonPath A character string specifying the path to the Python interpreter.
#' @param table Translation table, for which you can give an
#' NCBI genetic code number or name(https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
#' default 1.
#' @param exclude_stop Whether stop to extract amino acid when encounter stop codon.
#' c("yes","no"), default "yes".
#' @param ... Useless args.
#'
#' @return None.
#'
#' @export
setGeneric("fetch_sequence",
           function(object,
                    output_prefix = NULL,
                    pythonPath = NULL,
                    type = c("cds","exon"),
                    coding_type = c("NT","AA"),
                    table = 1,
                    exclude_stop = c("yes","no"), ...) standardGeneric("fetch_sequence"))




#' method for fetch_sequence
#'
#' @rdname fetch_sequence
#' @exportMethod fetch_sequence
setMethod("fetch_sequence",
          signature(object = "ribosomeObj"),
          function(object,
                   output_prefix = NULL,
                   pythonPath = NULL,
                   type = c("cds","exon"),
                   coding_type = c("NT","AA"),
                   table = 1,
                   exclude_stop = c("yes","no"),...){
            # check args
            type <- match.arg(type,c("cds","exon"))
            coding_type <- match.arg(coding_type,c("NT","AA"))
            exclude_stop <- match.arg(exclude_stop,c("yes","no"))

            # ============================================================================
            # run
            reticulate::py_config()
            if(reticulate::py_available() == FALSE){
              message("Please install python first!")
            }else{
              if(!is.null(pythonPath)){
                reticulate::use_python(pythonPath)
              }

              # check modual
              if (!reticulate::py_module_available("pyfaidx")) {
                reticulate::py_install("pyfaidx")
              }

              if (!reticulate::py_module_available("biopython")) {
                reticulate::py_install("biopython")
              }

              # run code
              pyscript.path = system.file("extdata", "fetchSequneceFromGeneAnno.py", package = "RiboProfiler")
              reticulate::source_python(pyscript.path)

              if(type == "cds" & coding_type == "NT"){
                output_file <- paste("annotation_data/",output_prefix,"_longest_cds.fa",sep = "")
                object@CDS.sequence <- output_file
              }

              if(type == "cds" & coding_type == "AA"){
                output_file <- paste("annotation_data/",output_prefix,"_longest_amino_acid.fa",sep = "")
                object@amino.acid.sequence <- output_file
              }

              suppressMessages(
                reticulate::py$fetchSequneceFromGeneAnno(gene_file = object@longest.anno.file,
                                                         genome_file = object@genome.file,
                                                         output_file = output_file,
                                                         type = type,
                                                         coding_type = coding_type,
                                                         table = table,
                                                         exclude_stop = exclude_stop)
              )

              return(object)
            }
          }
)
