#' Adjust Offset for RPF
#'
#' This function adjusts the offsets for quality control data based on the read lengths and associated offsets.
#' It separates the read lengths and offsets, performs necessary data type conversions, and then adjusts
#' the relative start positions by the absolute values of the offsets.
#'
#' @param offset_df A data frame containing the read lengths and offsets.
#' @param shift A shift for E/P/A site of reads position adjustment. Default 0.
#' @param object ribosomeObj object
#' @param ... useless args
#'
#' @return A data frame with adjusted relative start positions.
#'
#' @export
setGeneric("adjust_offset",function(object,offset_df = NULL,shift = 0,...) standardGeneric("adjust_offset"))





#' method for adjust_offset
#'
#' @rdname adjust_offset
#' @exportMethod adjust_offset
setMethod("adjust_offset",
          signature(object = "ribosomeObj"),
          function(object,
                   offset_df = NULL,
                   shift = 0,...){
            df_offset <- offset_df |>
              tidyr::separate_longer_delim(c(readLengths,Offsets),delim = ",") |>
              dplyr::rename(len = readLengths) |>
              dplyr::mutate(len = as.numeric(as.character(len)),
                            Offsets = as.numeric(Offsets))

            length_rpf <- unique(df_offset$len)

            shift_offset <- object@raw.counts |>
              dplyr::mutate(len = as.numeric(as.character(len))) |>
              dplyr::filter(len %in% length_rpf) |>
              dplyr::left_join(y = df_offset,by = c("sample","len")) |>
              dplyr::mutate(relst = relst + -(Offsets) + shift,
                            relsp = relsp + -(Offsets) + shift,
                            trans_pos = trans_pos + -(Offsets) + shift) |>
              dplyr::select(-Offsets,-bamFiles, -bamLegends)

            object@shifted.data <- shift_offset
            object@shifted <- "TRUE"

            return(object)
          }
)

