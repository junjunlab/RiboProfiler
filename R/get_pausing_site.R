#' Detect Pausing Sites
#'
#' This function processes normalized ribosome profiling data to detect pausing sites.
#' It writes intermediate normed data files and uses a Python script to perform
#' the actual detection of pausing sites. This method can be refered to pausepred
#' [PausePred and Rfeet: webtools for inferring ribosome pauses and visualizing footprint density from ribosome profiling data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140459/)
#'
#' @param object ribosomeObj object.
#' @param min_counts An integer specifying the minimum counts threshold for
#' considering a site as a pausing site. Default is \code{64}.
#' @param window An integer specifying the window size for detecting pausing
#' sites. Default is \code{1000}.
#' @param ... Useless args.
#'
#' @return This function does not return a value. It writes output files to the
#' "pausing_site_data" directory.
#'
#' @export
#' @examples
#' \dontrun{
#' # Example usage:
#' get_pausing_site(object = ribosomeObj,
#'                  min_counts = 64,
#'                  window = 500)
#' }
setGeneric("get_pausing_site",
           function(object,
                    min_counts = 64,
                    window = 1000, ...) standardGeneric("get_pausing_site"))





#' method for get_pausing_site
#'
#' @rdname get_pausing_site
#' @exportMethod get_pausing_site
setMethod("get_pausing_site",
          signature(object = "ribosomeObj"),
          function(object,
                   min_counts = 64,
                   window = 1000,...){
            # ========================================================================================
            # output normed data
            # ========================================================================================
            dir.create("pausing_site_data",showWarnings = F)

            normed_file <- object@normalized.data
            sp <- unique(normed_file$sample)

            lapply(seq_along(sp),function(x){
              tmp <- subset(normed_file,sample == sp[x])

              # output
              vroom::vroom_write(x = tmp,
                                 file = paste("pausing_site_data/",sp[x],"_normed.txt",sep = ""),
                                 col_names = F)

            }) -> none

            # ========================================================================================
            # detect pausing site
            # ========================================================================================
            # run code
            pyscript.path = system.file("extdata", "detectPausingSite.py", package = "RiboProfiler")
            reticulate::source_python(pyscript.path)

            # loop calculation
            purrr::map_df(seq_along(sp),function(x){
              normed_input_file <- paste("pausing_site_data/",sp[x],"_normed.txt",sep = "")
              output_file <- paste("pausing_site_data/",sp[x],"_pausingSite.txt",sep = "")

              # run
              suppressMessages(
                reticulate::py$detectPausingSite(longest_trans_file = object@longest.anno.file,
                                                 normed_file = normed_input_file,
                                                 output_file = output_file,
                                                 min_counts = as.integer(min_counts),
                                                 window = as.integer(window))
              )

              message(paste("#",normed_input_file,"has been processed!"))

              # ==========================================
              # load file
              ps <- vroom::vroom(file = output_file,col_names = FALSE,show_col_types = FALSE)
              colnames(ps) <- c("gene_name","trans_id","trans_pos","norm_exp","pausing_score")
              ps$sample <- sp[x]

              return(ps)
            }) -> df
          }
)
