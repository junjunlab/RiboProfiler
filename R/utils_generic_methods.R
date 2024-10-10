#' @title Subset data from ribosomeObj object.
#' @param object ribosomeObj object.
#' @param ... Useless args.
#'
#' @return ribosomeObj object.
#' @export
setGeneric("subdata", function(object, ...) standardGeneric("subdata"))



#' method for subdata
#'
#' @rdname subdata
#' @exportMethod subdata
setMethod("subdata",
          signature(object = "ribosomeObj"),
          function(object,...){
            raw_df <- object@raw.counts
            shift_df <- object@shifted.data
            norm_df <- object@normalized.data
            enrich_df <- object@enriched.data

            raw_df_ft <- base::subset(x = raw_df,...)
            object@raw.counts <- raw_df_ft

            if(nrow(shift_df) > 0){
              raw_shift_df <- base::subset(x = shift_df,...)
              object@shifted.data <- raw_shift_df
            }

            if(nrow(norm_df) > 0){
              raw_norm_df <- base::subset(x = norm_df,...)
              object@normalized.data <- raw_norm_df
            }

            if(nrow(enrich_df) > 0){
              raw_enrich_df <- base::subset(x = enrich_df,...)
              object@enriched.data <- raw_enrich_df
            }

            return(object)
          }
)





#' Generate Track Data Frame for Selected Genes
#'
#' This function extracts and prepares track data for selected genes from gene
#' annotation and normalized expression files.
#' It filters gene annotation data for selected genes, merges it with normalized
#' expression data, and formats it for further analysis or visualization.
#'
#' @param object ribosomeObj object.
#' @param select_gene A character vector of gene names to filter the gene annotation data.
#' If NULL, no filtering is applied. Default is NULL.
#' @param select_sample A character vector of sample names to filter. Default is NULL.
#' @param unit_to_codon Whether transform nuceotide position to codon position. Default is FALSE.
#' @param ... Useless args.
#'
#' @return A data frame containing the filtered and transformed data suitable for
#' tracking gene expression. The returned data frame includes columns for sample,
#' gene name, transcript ID, transcript position (transpos), normalized expression (density),
#' and a fixed column (type) set to 'ribo'.
#'
#' @importFrom dplyr filter select left_join rename mutate
#'
#' @examples
#' \dontrun{# Assuming `longest_trans_file` and `normed_file` are already defined:
#' track_df <- get_track_df(
#'   object = ribosomeObj,
#'   select_gene = c("Gene1", "Gene2")
#' )}
#'
#' @export
setGeneric("get_track_df",
           function(object,
                    select_gene = NULL,
                    select_sample = NULL,
                    unit_to_codon = FALSE, ...) standardGeneric("get_track_df"))




#' method for get_track_df
#'
#' @rdname get_track_df
#' @exportMethod get_track_df
setMethod("get_track_df",
          signature(object = "ribosomeObj"),
          function(object,
                   select_gene = NULL,
                   select_sample = NULL,
                   unit_to_codon = FALSE,...){

            # gene annotation
            gene_anao <- object@longest.annotation
            normed_file <- object@normalized.data

            colnames(gene_anao) <- c("id","gene_name","gene_id","trans_id","chrom","strand",
                                     "cds_rg","exon_rg","5UTR_length","CDS_length","3UTR_length")

            gene_selected <- gene_anao %>%
              dplyr::select(gene_name,trans_id,`5UTR_length`,CDS_length,`3UTR_length`) %>%
              dplyr::filter(gene_name %in% select_gene)

            # whether transform nucleotide to codon
            if(unit_to_codon == TRUE){
              mode = "codon"
              track_df <- normed_file %>%
                dplyr::filter(trans_id %in% gene_selected$trans_id) %>%
                dplyr::left_join(y = gene_selected,by = "trans_id") %>%
                dplyr::mutate(cdsft = dplyr::if_else(CDS_length%%3 == 0,1,0)) %>%
                # filter cds length %%3 == 0
                dplyr::filter(cdsft == 1) %>%
                dplyr::mutate(trans_pos = trans_pos - `5UTR_length`) %>%
                # transpos trans_pos into codon pos
                dplyr::mutate(codon_pos = dplyr::if_else(trans_pos%%3 == 0, trans_pos/3,trans_pos%/%3 + 1)) %>%
                dplyr::select(sample,gene_name,trans_id,codon_pos,norm_exp) %>%
                dplyr::rename(density = norm_exp,transpos = codon_pos) %>%
                dplyr::group_by(sample,gene_name,trans_id,transpos) %>%
                dplyr::summarise(density = sum(density)) %>%
                dplyr::mutate(type = "ribo")
            }else{
              mode = "nt"
              track_df <- normed_file %>%
                dplyr::filter(trans_id %in% gene_selected$trans_id) %>%
                dplyr::left_join(y = gene_selected,by = "trans_id") %>%
                dplyr::select(sample,gene_name,trans_id,trans_pos,norm_exp) %>%
                dplyr::rename(density = norm_exp,transpos = trans_pos) %>%
                dplyr::group_by(sample,gene_name,trans_id,transpos) %>%
                dplyr::summarise(density = sum(density)) %>%
                dplyr::mutate(type = "ribo")
            }

            if(!is.null(select_sample)){
              track_df <- subset(track_df,sample %in% select_sample)
            }else{
              track_df
            }

            # assign a unit attribute
            attr(track_df,"unit") <- mode

          }
)
