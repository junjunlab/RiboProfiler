#' Adjust Offset for RPF
#'
#' This function adjusts the offsets for quality control data based on the read lengths and associated offsets.
#' It separates the read lengths and offsets, performs necessary data type conversions, and then adjusts
#' the relative start positions by the absolute values of the offsets.
#'
#' @param offset_df A data frame containing the read lengths and offsets.
#' @param qc_df A data frame containing quality control data with sample, length,
#' relative start positions, and other relevant information.
#' @param shift A shift for E/P/A site of reads position adjustment. Default 0.
#'
#' @return A data frame with adjusted relative start positions.
#'
#' @export
adjust_offset <- function(offset_df = NULL,
                          qc_df = NULL,
                          shift = 0){
  df_offset <- offset_df |>
    tidyr::separate_longer_delim(c(readLengths,Offsets),delim = ",") |>
    dplyr::rename(length = readLengths) |>
    dplyr::mutate(length = as.numeric(length),
                  Offsets = abs(as.numeric(Offsets)) + shift)

  length_rpf <- unique(df_offset$length)

  shift_offset <- qc_df |>
    dplyr::mutate(length = as.numeric(as.character(length))) |>
    dplyr::filter(length %in% length_rpf) |>
    dplyr::left_join(y = df_offset,by = c("sample","length")) |>
    dplyr::mutate(relst = relst + abs(Offsets),
                  relsp = relsp + abs(Offsets),
                  trans_pos = trans_pos + abs(Offsets)) |>
    dplyr::select(-Offsets,-bamFiles, -bamLegends)

  return(shift_offset)
}
