globalVariables(c("bam_fie", "bam_file", "maxpos", "relst"))

#' psite_offset_check
#'
#' @param object ribosomeObj object.
#' @param relative_distance A vector of two values representing the range of relative distances
#'   to filter by. Default is c(-18,12).
#' @param read_length A vector of two values representing the range of read lengths to filter by.
#'   Default is c(27,32).
#' @param label_size The size of the text labels. Default is 3.
#' @param label_pos A vector of two values representing the position of the labels. Default is c(0.9,0.9).
#' @param line_col The color of the line. Default is "#993366".
#' @param peak_line_col The color of the peak line. Default is "red".
#' @param ... Useless args.
#'
#' @return A list with two elements:
#' \item{summary_offset}{A data frame containing the summary of offsets.}
#' \item{plot}{A ggplot object representing the plot of read counts against distance.}
#'
#' @export
setGeneric("psite_offset_check",
           function(object,
                    relative_distance = c(-20,20),
                    read_length = c(27,32),
                    label_size = 3,
                    label_pos = c(0.9,0.9),
                    line_col = "#993366",
                    peak_line_col = "red",...) standardGeneric("psite_offset_check"))




#' method for psite_offset_check
#'
#' @rdname psite_offset_check
#' @exportMethod psite_offset_check
setMethod("psite_offset_check",
          signature(object = "ribosomeObj"),
          function(object,
                   relative_distance = c(-20,20),
                   read_length = c(27,32),
                   label_size = 3,
                   label_pos = c(0.9,0.9),
                   line_col = "#993366",
                   peak_line_col = "red",...){
            # ===================================================================================
            # process data
            # ===================================================================================
            qc_data <- object@raw.counts
            qc_data$length <- as.numeric(as.character(qc_data$length))

            df_ft <- qc_data |>
              dplyr::filter(relst >= relative_distance[1] & relst <= relative_distance[2]) |>
              dplyr::filter(length >= read_length[1] & length <= read_length[2]) |>
              dplyr::group_by(sample,length,relst) |>
              dplyr::summarise(count = sum(counts))

            # get max peak position
            df_maxht <- df_ft |>
              dplyr::group_by(sample,length) |>
              dplyr::summarise(maxpos = max(count)) |>
              dplyr::left_join(y = df_ft,by = c("sample","length"),
                               multiple = "all") |>
              dplyr::filter(maxpos == count)

            # retain only one offset
            purrr::map_df(unique(df_maxht$sample),function(x){
              tmp <- subset(df_maxht,sample == x)
              len <- unique(tmp$length)
              purrr::map_df(seq_along(len),function(x){
                tmp2 <- subset(tmp,length == len[x]) %>%
                  dplyr::arrange(relst) %>%
                  head(n = 1)
                return(tmp2)
              }) -> uniq_offset
              return(uniq_offset)
            }) ->df_maxht

            # summarise offsets
            summary_offset <-
              plyr::ddply(df_maxht[,c("sample","length","relst")],
                          plyr::.(sample), dplyr::mutate,
                          readLengths = paste(length, collapse = ","),
                          Offsets = paste(relst, collapse = ",")) |>
              dplyr::select(-length,-relst)|> unique() |>
              dplyr::mutate(bamFiles = sample,bamLegends = sample)

            # ===================================================================================
            # plot data
            # ===================================================================================
            p <-
              ggplot(df_ft) +
              geom_line(aes(x = relst,y = count),linewidth = 0.75,color = line_col) +
              geom_vline(xintercept = 0,
                         lty = "dashed",color = peak_line_col,linewidth = 0.75) +
              geom_vline(data = df_maxht,
                         aes(xintercept = relst),
                         lty = "dashed",color = peak_line_col,linewidth = 0.75) +
              ggpp::geom_text_npc(data = df_maxht,
                                  aes(npcx = label_pos[1],npcy = label_pos[2],
                                      label = paste(abs(relst),"nt",sep = " ")),
                                  size = label_size,
                                  fontface = "bold") +
              theme_bw() +
              theme(strip.text.y = element_text(angle = 0),
                    axis.text = element_text(colour = "black"),
                    panel.grid = element_blank(),
                    strip.text = element_text(face = "bold.italic",size = rel(1))) +
              # facet_grid(length~sample,scales = "free_y") +
              ggh4x::facet_grid2(length~sample,scales = "free",independent = "y") +
              ylab("Read counts") +
              xlab("Distance of 5'end to start codon (nt)")

            # ====================================================================================
            # return
            return(list(summary_offset = summary_offset,
                        plot = p))
          }
)
