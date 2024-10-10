#' qc_plot function
#'
#' A function to create quality control (QC) plots for read length, frame and
#' features.
#'
#' @param object ribosomeObj object.
#' @param assay Using raw data or shifted data to plot, c("raw","shifted").
#' @param type The type of QC plot to create. Must be one of "length",
#' "length_frame","frame" or "feature".
#' @param geom_col_list A list of additional ggplot2 layers to add to the plot
#' using the \code{geom_col} function.
#' @param facet_wrap_list A list of additional ggplot2 layers to add to the plot
#' using the \code{facet_wrap} function.
#' @param ... Useless args.
#'
#' @return Returns a ggplot2 object representing the created QC plot.
#'
#' @examples
#' \dontrun{
#' qc_plot(qc_data = my_data, type = "length")
#' }
#'
#' @import ggplot2
#' @importFrom dplyr filter mutate select group_by summarise
#' @import utils
#'
#' @export
setGeneric("qc_plot",
           function(object,
                    assay = c("raw","shifted"),
                    type = c("length","length_frame","feature","frame"),
                    geom_col_list = list(),
                    facet_wrap_list = list(),...) standardGeneric("qc_plot"))





#' method for qc_plot
#'
#' @rdname qc_plot
#' @exportMethod qc_plot
setMethod("qc_plot",
          signature(object = "ribosomeObj"),
          function(object,
                   assay = c("raw","shifted"),
                   type = c("length","length_frame","feature","frame"),
                   geom_col_list = list(),
                   facet_wrap_list = list(),...){
            assay <- match.arg(assay,c("raw","shifted"))
            options(warn=-1)
            # Suppress summarise info
            options(dplyr.summarise.inform = FALSE)

            # ===================================================================
            # plot layers
            # ===================================================================
            # plot
            pmain <-
              ggplot() +
              theme_bw() +
              jj_theme() +
              xlab('Read length') + ylab('Reads numbers(k)') +
              do.call(facet_wrap,
                      modifyList(list(~sample,scales = 'free',ncol = 2),
                                 facet_wrap_list))

            # process layers
            if(assay == "raw"){
              qc_data <- object@raw.counts
            }else{
              qc_data <- object@shifted.data
            }

            type <- match.arg(type,c("len","length_frame","feature","frame"))
            if(type == "length"){
              len <- qc_data %>% group_by(group,sample,len) %>%
                dplyr::summarise(num = sum(counts))

              layer_tmp <- pmain +
                do.call(geom_col,
                        modifyList(list(data = len,
                                        mapping = aes(x = len,y = num/1000),
                                        fill = "#A4BE7B",width = 0.6),
                                   geom_col_list))
            }else if(type == "length_frame"){
              frame <- qc_data %>% group_by(group,sample,len,framest) %>%
                dplyr::summarise(num = sum(counts))

              layer_tmp <- pmain +
                do.call(geom_col,
                        modifyList(list(data = frame,
                                        mapping = aes(x = len,y = num/1000,fill = factor(framest)),
                                        position = position_dodge2()),
                                   geom_col_list)) +
                scale_fill_brewer(palette = 'Greens',name = '',direction = -1,
                                  labels = c('frame 0','frame 1','frame 2'))
            }else if(type == "feature"){
              featuredf <- qc_data %>% group_by(group,sample,feature) %>%
                dplyr::summarise(num = sum(counts))

              layer_tmp <- pmain +
                do.call(geom_col,
                        modifyList(list(data = featuredf,
                                        mapping = aes(x = rev(factor(feature)),y = num/1000,fill = factor(feature)),
                                        width = 0.6,
                                        position = position_dodge2()),
                                   geom_col_list)) +
                scale_x_discrete(labels = c('CDS',"5'UTR","3'UTR")) +
                scale_fill_brewer(palette = 'Blues',name = '',direction = 1,
                                  labels = c("3'UTR","5'UTR",'CDS'))
            }else if(type == "frame"){
              frame <- qc_data %>% group_by(group,sample,framest) %>%
                dplyr::summarise(num = sum(counts))

              layer_tmp <- pmain +
                do.call(geom_col,
                        modifyList(list(data = frame,
                                        mapping = aes(x = factor(framest),y = num/1000,fill = factor(framest)),
                                        width = 0.6,
                                        position = position_dodge2()),
                                   geom_col_list)) +
                scale_fill_brewer(palette = 'Greens',name = '',direction = -1,
                                  labels = c('frame 0','frame 1','frame 2'))
            }

            return(layer_tmp)
          }
)







#' Plot the distribution of reads relative to start or stop codon
#'
#' This function takes in a data frame containing quality control data for a set
#' of samples, and plots the
#' distribution of reads relative to either the start or stop codon. The user can
#' specify the type of distribution,as well as a distance range to filter the data.
#' The function returns a ggplot object.
#'
#' @param object ribosomeObj object.
#' @param assay Using raw data or shifted data to plot, c("raw","shifted").
#' @param type A character string specifying the type of distribution. Must be
#' one of "relst" (relative to start codon) or "relsp" (relative to stop codon).
#' @param dist_range A numeric vector specifying the minimum and maximum distances
#' to include in the plot. By default, the entire range(-50-100/-100-50) is included.
#' @param shift A numeric value specifying the p site to shift, default 0.
#' @param geom_col_list A list of arguments to pass to the
#' \code{\link[ggplot2]{geom_col}} layer.
#' @param facet_wrap_list A list of arguments to pass to the
#' \code{\link[ggplot2]{facet_wrap}} layer.
#' @param ... Useless args.
#'
#' @importFrom rlang ensyms
#'
#' @return A ggplot object displaying the distribution of reads relative to the
#' start or stop codon.
#'
#' @examples
#' \dontrun{
#' rel_to_start_stop(qc_data = my_qc_data, type = "relst", dist_range = c(-50, 100))
#' }
#' @export
setGeneric("rel_to_start_stop",
           function(object,
                    assay = c("raw","shifted"),
                    type = c("relst","relsp"),
                    dist_range = NULL,
                    shift = 0,
                    geom_col_list = list(),
                    facet_wrap_list = list(),...) standardGeneric("rel_to_start_stop"))







#' method for rel_to_start_stop
#'
#' @rdname rel_to_start_stop
#' @exportMethod rel_to_start_stop
setMethod("rel_to_start_stop",
          signature(object = "ribosomeObj"),
          function(object,
                   assay = c("raw","shifted"),
                   type = c("relst","relsp"),
                   dist_range = NULL,
                   shift = 0,
                   geom_col_list = list(),
                   facet_wrap_list = list(),...){
            assay <- match.arg(assay,c("raw","shifted"))
            options(warn=-1)
            # Suppress summarise info
            options(dplyr.summarise.inform = FALSE)

            # ============================================================================
            # process arguments
            # ============================================================================
            type <- match.arg(type,c("relst","relsp"))

            if(type == "relst"){
              v = "framest"
              vars_f <- rlang::ensyms(type,v)
              d_range <- c(-50,100)
            }else{
              v = "framesp"
              vars_f <- rlang::ensyms(type,v)
              d_range <- c(-100,50)
            }

            if(!is.null(dist_range)){
              d_range <- dist_range
            }
            # ============================================================================
            # process data
            # ============================================================================
            if(assay == "raw"){
              qc_data <- object@raw.counts
            }else{
              qc_data <- object@shifted.data
            }


            rel <- qc_data %>% group_by(sample,!!!vars_f) %>%
              dplyr::summarise(allCounts = sum(counts))

            # total reads
            totalReads <- qc_data %>% group_by(sample) %>%
              dplyr::summarise(total = sum(counts))

            # normalize to total reads
            s = unique(rel$sample)
            plyr::ldply(s,function(x){
              tmp <- rel %>% filter(sample == x)
              total_ct <- unlist(totalReads[which(totalReads$sample == x),"total"]) %>% as.numeric()
              tmp$normaledReads <- (tmp$allCounts/total_ct)*10^6
              return(tmp)
            }) -> rel

            # filter region
            df <- rel %>% filter(!!vars_f[[1]] >= d_range[1] & !!vars_f[[1]] <= d_range[2])
            df <- df %>% mutate(!!vars_f[[2]]  := factor(!!vars_f[[2]] ,levels = c(0,1,2)))

            # ============================================================================
            # plot
            # ============================================================================
            pmain <-
              ggplot() +
              theme_bw() +
              jj_theme() +
              xlab('Read length') + ylab('Reads numbers(k)') +
              do.call(facet_wrap,
                      modifyList(list(~sample,scales = 'free',ncol = 2),
                                 facet_wrap_list)) +
              do.call(geom_col,modifyList(
                list(data = df,
                     mapping = aes(x = !!vars_f[[1]] + shift,y = normaledReads/1000,
                                   fill = factor(!!vars_f[[2]])),
                     width = 0.8,
                     position = position_dodge2()),
                geom_col_list)) +
              scale_fill_brewer(palette = 'OrRd',name = '',direction = -1,
                                labels = c('frame 0','frame 1','frame 2'))

            # layers
            if(type == "relst"){
              pfin <-
                pmain +
                xlab('Relative to start codon') +
                ylab('Normalizee Reads counts(k)')
            }else{
              pfin <-
                pmain +
                xlab('Relative to stop codon') +
                ylab('Normalizee Reads counts(k)')
            }

            return(pfin)
          }
)
