globalVariables(c("Amplitude", "Period"))

#' Check Periodicity in Ribosome Profiling Data
#'
#' This function performs a periodicity check on reads from ribosome profiling data.
#' It uses the Discrete Fourier Transform to evaluate the periodicity of ribosomal reads
#' around start codons, and visualizes the result.
#'
#' @param object A `ribosomeObj` object containing the raw counts and metadata.
#' @param read_length A numeric vector of length 2 specifying the range of read lengths to consider.
#'   Default is c(27, 30).
#' @param relative_dist A numeric vector of length 2 indicating the distance range around the start codon to evaluate.
#'   Default is c(0, 150).
#' @param period_max A numeric value indicating the maximum period (in nucleotides) to display in the plot.
#'   Default is 10.
#' @param merge_rep Whether merge replicates. Default is FALSE.
#' @param ... Additional arguments to be passed to or from other methods.
#'
#' @return A ggplot object displaying the amplitude of different periods across samples and read lengths.
#'
#' @importFrom stats fft
#' @importFrom scales label_math
#' @importFrom ggh4x facet_grid2
#'
#'
#' @export
setGeneric("periodicity_check",
           function(object,
                    read_length = c(27,30),
                    relative_dist = c(0,150),
                    period_max = 10,
                    merge_rep = FALSE,...) standardGeneric("periodicity_check"))




#' method for periodicity_check
#'
#' @rdname periodicity_check
#' @exportMethod periodicity_check
setMethod("periodicity_check",
          signature(object = "ribosomeObj"),
          function(object,
                   read_length = c(27,30),
                   relative_dist = c(0,150),
                   period_max = 10,
                   merge_rep = FALSE,
                   ...){

            raw <- object@raw.counts

            sp <- unique(raw$sample)

            # loop for sample
            # x = 1
            purrr::map_df(seq_along(sp),function(x){
              tmp <- subset(raw,sample == sp[x])

              anno_raw <- tmp %>%
                dplyr::filter(len %in% c(read_length[1]:read_length[2])) %>%
                dplyr::group_by(relst,len) %>%
                dplyr::summarise(nc = sum(counts)) %>%
                dplyr::filter(relst >= relative_dist[1] & relst <= relative_dist[2])

              # loop for length
              len <- unique(anno_raw$len)
              purrr::map_df(seq_along(len),function(l){
                tmp2 <- subset(anno_raw,len == len[l])

                # Discrete Fourier transform of reads
                freq_domain_data <- fft(tmp2$nc)
                power_spectrum <- Mod(freq_domain_data)^2
                N <- length(freq_domain_data)

                periods <- N / seq(1, N/2)

                # to data frame
                df <- data.frame(Period = periods,
                                 Amplitude = power_spectrum[2:(N/2 + 1)],
                                 length = len[l],
                                 rep = tmp$rep[1],
                                 sample = sp[x]) %>%
                  dplyr::filter(Period <= period_max)

                return(df)
              }) -> tmp_df
              return(tmp_df)
            }) -> fft_res

            # ==================================================================
            # merge replicates
            if(merge_rep == TRUE){
              fft_res <- fft_res %>%
                dplyr::group_by(Period,length,rep) %>%
                dplyr::summarise(Amplitude = mean(Amplitude)) %>%
                dplyr::rename(sample = rep)
            }

            # ==================================================================
            # plot

            ggplot(fft_res) +
              geom_line(aes(x = Period,y = Amplitude,color = sample),
                        linewidth = 0.5) +
              theme_bw() +
              # facet_grid(sample~length,scales = "free",axes = "all") +
              ggh4x::facet_grid2(sample~length,scales = "free",axes = "all",independent = "all") +
              theme(axis.text = element_text(colour = "black"),
                    panel.grid = element_blank(),
                    strip.text = element_text(face = "bold",size = rel(1)),
                    strip.text.x = element_text(face = "bold.italic",size = rel(1)),
                    strip.text.y.right = element_text(angle = 0,hjust = 0),
                    strip.background = element_blank()) +
              scale_y_continuous(labels = scales::label_math()) +
              xlab("Period | around start codon (nt)")
          })
