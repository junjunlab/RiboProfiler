globalVariables(c("counts", "dfqc", "feature", "framest", "group", "num",
                  "normaledReads",":=","density", "gene_name", "if_else",
                  "max_v", "min_v", "remove_chrom_panel_border", "rg_label",
                  "trans_id", "transpos", "type", "ymax", "ymin",".","region",
                  "id","calculateRibosomeDensity","map_type", "reads","end",
                  "gene_id", "seqnames", "start", "strand", "transcript_id",
                  "utr3"))

#' theme setting
#' Theme for plots
#'
#' @return theme list
#' @export
jj_theme <- function(){
  theme(strip.background = element_rect(),
        axis.text = element_text(color = "black",size = rel(1)),
        axis.title = element_text(color = "black",size = rel(1.2),face = "bold"),
        strip.text = element_text(color = "black",size = rel(1),face = "bold"),
        strip.text.x = element_text(color = "black",size = rel(1),face = "bold.italic"),
        panel.grid = element_blank())
}

#' qc_plot function
#'
#' A function to create quality control (QC) plots for read length, frame and
#' features.
#'
#' @param qc_data The input data frame for QC analysis. Requires columns "group",
#' "sample", "length", "counts", "framest" and "feature".
#' @param type The type of QC plot to create. Must be one of "length",
#' "length_frame","frame" or "feature".
#' @param geom_col_list A list of additional ggplot2 layers to add to the plot
#' using the \code{geom_col} function.
#' @param facet_wrap_list A list of additional ggplot2 layers to add to the plot
#' using the \code{facet_wrap} function.
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
qc_plot <- function(qc_data = NULL,
                    type = c("length","length_frame","feature","frame"),
                    geom_col_list = list(),
                    facet_wrap_list = list()){
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
  type <- match.arg(type,c("length","length_frame","feature","frame"))
  if(type == "length"){
    len <- qc_data %>% group_by(group,sample,length) %>%
      dplyr::summarise(num = sum(counts))

    layer_tmp <- pmain +
      do.call(geom_col,
              modifyList(list(data = len,
                              mapping = aes(x = length,y = num/1000),
                              fill = "#A4BE7B",width = 0.6),
                         geom_col_list))
  }else if(type == "length_frame"){
    frame <- qc_data %>% group_by(group,sample,length,framest) %>%
      dplyr::summarise(num = sum(counts))

    layer_tmp <- pmain +
      do.call(geom_col,
              modifyList(list(data = frame,
                              mapping = aes(x = length,y = num/1000,fill = factor(framest)),
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


#' Plot the distribution of reads relative to start or stop codon
#'
#' This function takes in a data frame containing quality control data for a set
#' of samples, and plots the
#' distribution of reads relative to either the start or stop codon. The user can
#' specify the type of distribution,as well as a distance range to filter the data.
#' The function returns a ggplot object.
#'
#' @param qc_data A data frame containing quality control data for a set of samples.
#' @param type A character string specifying the type of distribution. Must be
#' one of "relst" (relative to start codon) or "relsp" (relative to stop codon).
#' @param dist_range A numeric vector specifying the minimum and maximum distances
#' to include in the plot. By default, the entire range(-50-100/-100-50) is included.
#' @param shift A numeric value specifying the p site to shift, default 0.
#' @param geom_col_list A list of arguments to pass to the
#' \code{\link[ggplot2]{geom_col}} layer.
#' @param facet_wrap_list A list of arguments to pass to the
#' \code{\link[ggplot2]{facet_wrap}} layer.
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
rel_to_start_stop <- function(qc_data = NULL,
                              type = c("relst","relsp"),
                              dist_range = NULL,
                              shift = 0,
                              geom_col_list = list(),
                              facet_wrap_list = list()){
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



#' Plot ribosome density and RNA coverage along transcript
#'
#' This function plots ribosome density and RNA coverage along the transcript,
#' with gene structures overlaid. The input data should be a data frame with
#' columns "transpos", "density", "type" (either "ribo" or "rna"), "gene_name",
#' "sample". Gene structures are loaded from a file named "longest_info.txt".
#' You can customize various aspects of the plot, such as color scheme, panel
#' size, and panel borders.
#'
#' @param signal_data A data frame with columns "transpos", "density", "type"
#' (either "ribo" or "rna"), "gene_name", "sample".
#' @param gene_anno Longest transcript gene annotation file.
#' @param plot_type the plot type for track, "translatome"(default) or "interactome".
#' @param add_gene_struc Whether add gene structure plot. Default is TRUE.
#' @param show_cds_region_only Whether retain CDS region only. Default is FALSE.
#' @param structure_col Color of the gene structure rectangles.
#' @param background_col Color of the plot background. Default is "grey90".
#' @param range_pos Position of the range labels. Default is c(0.85,0.85).
#' @param reverse_rna Whether reverse RNA track. Default is TRUE.
#' @param rna_signal_scale The scale value of rna coverage track. Default is 1.
#' It is used for better visualization when the difference of ribo density value
#' and rna coverage value is quite huge.
#' @param signal_col A named vector of colors for ribosome and RNA signals.
#' Default is NULL, which means using a default color scheme.
#' @param line_col A named vector of colors for enrichment line signals.
#' Default is NULL, which means using a default color scheme.
#' @param panel_size A numeric vector of length 2 specifying the width and height
#' of each panel. Default is NULL, which means using a default panel size.
#' @param remove_all_panel_border Logical indicating whether to remove all panel
#' borders. Default is FALSE.
#' @param remove_trans_panel_border Logical indicating whether to remove the
#' border of the transcript panel only. Default is TRUE.
#' @param gene_order An optional character vector specifying the order of genes
#' in the plot. By default, genes are ordered by their appearance in the input data.
#' @param sample_order An optional character vector specifying the order of
#' samples in the plot. By default, samples are ordered by their appearance in
#' the input data.
#' @param fixed_col_range Logical indicating whether to use a fixed color range
#' for all panels. Default is TRUE.
#' @param show_ribo_only Show only ribo density track. Default is FALSE.
#' @param show_x_ticks show_x_ticks Whether show x aixs ticks. Default is FALSE.
#' @param show_y_ticks Whether show y aixs ticks. Default is FALSE.
#' @param sample_group_info The groups for samples, giving a named list with
#' samples, default NULL.
#' @param geom_col_params The parameters passed by ggplot2::geom_col function for
#' "translatome" plot type to overwrite params.
#' @param geom_line_params The parameters passed by ggplot2::geom_line function for
#' "interactome" plot type to overwrite params.
#' @param geom_hline_params The parameters passed by ggplot2::geom_hline function for
#' "interactome" plot type to overwrite params.
#' @param geom_ribbon_params The parameters passed by ggplot2::geom_ribbon function for
#' "interactome" plot type to overwrite params.
#' @param structure_label_y The y position of gene structure, defalut -3.
#' @param structure_label_size The text size of gene structure labels, defalut 2.5.
#' @param signal_range_digits The signal range digits, you can set it bigger when
#' signal range is small than 1, defalut 0.
#'
#' @importFrom zplyr geom_abs_text
#' @importFrom ggh4x facet_nested facetted_pos_scales force_panelsizes
#' @importFrom grid grid.ls grid.draw grid.newpage
#'
#' @export
track_plot <- function(signal_data = NULL,
                       gene_anno = NULL,
                       plot_type = c("translatome","interactome"),
                       add_gene_struc = TRUE,
                       show_cds_region_only = FALSE,
                       structure_label_y = -3,
                       structure_label_size = 2.5,
                       structure_col = NULL,
                       background_col = "white",
                       range_pos = c(0.8,0.85),
                       signal_range_digits = 0,
                       reverse_rna = TRUE,
                       rna_signal_scale = 1,
                       signal_col = NULL,
                       line_col = NULL,
                       panel_size = NULL,
                       remove_all_panel_border = FALSE,
                       remove_trans_panel_border = TRUE,
                       gene_order = NULL,
                       sample_order = NULL,
                       fixed_col_range = TRUE,
                       show_ribo_only = FALSE,
                       show_x_ticks = FALSE,
                       show_y_ticks = FALSE,
                       sample_group_info = NULL,
                       geom_col_params = list(),
                       geom_line_params = list(),
                       geom_hline_params = list(),
                       geom_ribbon_params = list()){
  # ==============================================================================
  # process data
  # ==============================================================================
  # show_ribo_only = T
  if(show_ribo_only == TRUE){
    signal_data <- signal_data %>% dplyr::filter(type != "rna")
  }else{
    if(reverse_rna == TRUE){
      signal_data <- signal_data %>%
        dplyr::mutate(density = dplyr::if_else(type %in% "rna",
                                               -density*rna_signal_scale,
                                               density))
    }else{
      signal_data <- signal_data %>%
        dplyr::mutate(density = dplyr::if_else(type %in% "rna",
                                               density*rna_signal_scale,
                                               density))
    }

    # order
    signal_data$type <- factor(signal_data$type,levels = c('rna','ribo'))
  }

  # add group info for samples
  # sample_group_info = list(group1 = c("sample1","sample2"))
  if(!is.null(sample_group_info)){
    # x = 1
    plyr::ldply(1:length(sample_group_info),function(x){
      tmp <- signal_data %>%
        dplyr::filter(sample %in% sample_group_info[[x]]) %>%
        dplyr::mutate(group = names(sample_group_info)[x])
      return(tmp)
    }) -> signal_data
  }

  # ==============================================================================
  # change trans_pos to cds_pos
  # ==============================================================================
  if(show_cds_region_only == TRUE){
    # load geneinfo
    geneInfo <- read.table(gene_anno)
    colnames(geneInfo) <- c('id','gene_name','gene_id','trans_id','chr','strand',
                            'cds_region','exon_region','utr5','cds','utr3')
    geneInfo <- geneInfo %>%
      dplyr::filter(gene_name %in% unique(signal_data$gene_name)) %>%
      dplyr::select(gene_name,utr5,cds,utr3)

    signal_data <- signal_data %>%
      dplyr::left_join(y = geneInfo,by = "gene_name") %>%
      dplyr::mutate(transpos = transpos - utr5) %>%
      dplyr::filter(transpos > 0 & transpos <= cds)
  }

  # ============================================================================
  # sample or gene orders
  # ============================================================================
  # gene_order = NULL
  # sample_order = NULL

  if(is.null(gene_order)){
    signal_data$gene_name <- factor(signal_data$gene_name,
                                    levels = unique(signal_data$gene_name))
  }else{
    signal_data$gene_name <- factor(signal_data$gene_name,
                                    levels = gene_order)
  }

  if(is.null(sample_order)){
    signal_data$sample <- factor(signal_data$sample,
                                 levels = unique(signal_data$sample))
  }else{
    signal_data$sample <- factor(signal_data$sample,
                                 levels = sample_order)
  }


  # ==============================================================================
  # gene strctures
  # ==============================================================================
  plot_type <- match.arg(plot_type,c("translatome","interactome"))

  if(plot_type == "interactome"){
    gene <- unique(signal_data$gene_name)

    # x = 1
    structure_df <- plyr::ldply(seq_along(gene), function(x){
      tmp <- signal_data[which(signal_data$gene_name %in% gene[x]),]

      df <- data.frame(gene_name = gene[x],
                       start = 1,
                       end = length(unique(tmp$codon_pos)),
                       ymin = -1,
                       ymax = 1,
                       sample = "trans",
                       region = "CDS",
                       group = NA)

      return(df)
    })
  }else if(plot_type == "translatome"){
    # load geneinfo
    geneInfo <- read.table(gene_anno)
    colnames(geneInfo) <- c('id','gene_name','gene_id','trans_id','chr','strand',
                            'cds_region','exon_region','utr5','cds','utr3')
    geneInfo <- geneInfo %>%
      dplyr::filter(gene_name %in% unique(signal_data$gene_name))

    # x = 1
    structure_df <- plyr::ldply(1:nrow(geneInfo), function(x){
      tmp <- geneInfo[x,]

      if(show_cds_region_only == TRUE){
        df <- data.frame(gene_name = tmp$gene_name,
                         start = 1,
                         end = tmp$cds,
                         ymin = -1,
                         ymax = 1,
                         sample = "trans",
                         region = "CDS",
                         group = NA)
      }else{
        df <- data.frame(gene_name = tmp$gene_name,
                         start = c(0,tmp$utr5,tmp$utr5 + tmp$cds),
                         end = c(tmp$utr5,tmp$utr5 + tmp$cds,tmp$utr5 + tmp$cds + tmp$utr3),
                         ymin = c(-0.5,-1,-0.5),
                         ymax = c(0.5,1,0.5),
                         sample = "trans",
                         region = c("5UTR","CDS","3UTR"),
                         group = NA)
      }

      return(df)
    })
  }

  # reassign orders
  structure_df$sample <- factor(structure_df$sample ,levels = c(levels(signal_data$sample),"trans"))
  structure_df$gene_name <- factor(structure_df$gene_name ,levels = levels(signal_data$gene_name))

  # ==============================================================================
  # panel range settings
  # ==============================================================================
  n_sample = length(unique(unique(signal_data$sample)))
  n_gene = length(unique(unique(signal_data$gene_name)))

  # fixed_col_range = F
  if(fixed_col_range == TRUE){
    track_range <- signal_data %>%
      dplyr::group_by(gene_name)

    # check plot_type
    if(plot_type == "interactome"){
      track_range <- track_range %>%
        dplyr::summarise(min_v = round(min(density - density_sd),digits = signal_range_digits),
                         max_v = round(max(density + density_sd),digits = signal_range_digits))
    }else if(plot_type == "translatome"){
      track_range <- track_range %>%
        dplyr::summarise(min_v = round(min(density),digits = signal_range_digits),
                         max_v = round(max(density),digits = signal_range_digits))
    }


    # reassign minimum value
    if(show_ribo_only == TRUE | reverse_rna == FALSE | plot_type == "interactome"){
      track_range$min_v <- 0
    }

    track_range <- track_range %>%
      dplyr::mutate(rg_label = paste("[",min_v,"-",max_v,"]",sep = "")) %>%
      replicate(n_sample,.,simplify = F) %>%
      do.call("rbind",.)

    track_range$sample <- rep(unique(signal_data$sample),each = n_gene)
  }else{
    track_range <- signal_data %>%
      dplyr::group_by(gene_name,sample)

    # check plot_type
    if(plot_type == "interactome"){
      track_range <- track_range %>%
        dplyr::summarise(min_v = round(min(density - density_sd),digits = 2),
                         max_v = round(max(density + density_sd),digits = 2))
    }else if(plot_type == "translatome"){
      track_range <- track_range %>%
        dplyr::summarise(min_v = round(min(density),digits = 2),
                         max_v = round(max(density),digits = 2))
    }

    # reassign minimum value
    if(show_ribo_only == TRUE | reverse_rna == FALSE | plot_type == "interactome"){
      track_range$min_v <- 0
    }

    track_range <- track_range %>%
      dplyr::mutate(rg_label = paste("[",min_v,"-",max_v,"]",sep = ""))
  }

  # add group info
  if(!is.null(sample_group_info)){
    # x = 1
    plyr::ldply(1:length(sample_group_info),function(x){
      tmp <- track_range %>%
        dplyr::filter(sample %in% sample_group_info[[x]]) %>%
        dplyr::mutate(group = names(sample_group_info)[x])
      return(tmp)
    }) -> track_range
  }

  track_range$sample <- factor(track_range$sample ,levels = levels(signal_data$sample))
  track_range$gene_name <- factor(track_range$gene_name ,levels = levels(signal_data$gene_name))

  # orders
  track_range <- track_range[order(track_range$sample,track_range$gene_name),]

  # y axis range layer
  range_layer <- lapply(1:(nrow(track_range) + length(unique(track_range$gene_name))), function(x){
    position <- ifelse(show_y_ticks == TRUE,"right","left")
    if(x <= nrow(track_range)){
      tmp <- track_range[x,]
      scale_y_continuous(limits = c(tmp$min_v,tmp$max_v),position = position)
    }else{
      scale_y_continuous(limits = c(-3,1),breaks = NULL,position = position)
    }
  })

  # ==============================================================================
  # paras settings
  # ==============================================================================

  # color
  if(is.null(signal_col)){
    signal_col <- c("ribo" = "#D21312","rna" = "#009FBD")
  }else{
    signal_col <- signal_col
  }

  # interactome line color
  if(is.null(line_col)){
    line_col_layer <- list(scale_color_brewer(palette = "Paired",name = ""),
                           scale_fill_brewer(palette = "Paired",name = ""))

  }else{
    line_col_layer <- list(scale_color_manual(values = line_col,name = ""),
                           scale_fill_manual(values = line_col,name = ""))
  }

  # structure_col = NULL
  if(is.null(structure_col)){
    structure_col <- c("5UTR" = "grey70","CDS" = "grey50","3UTR" = "grey70")
  }else{
    structure_col <- structure_col
  }

  # panel size
  if(is.null(panel_size)){
    panel_size <- c(4,1)
  }else{
    panel_size <- panel_size
  }

  # facet label
  label_df <- signal_data %>% dplyr::ungroup() %>%
    dplyr::select(gene_name,trans_id) %>% unique()

  new_label <- paste(label_df$gene_name,label_df$trans_id,sep = "\n")
  names(new_label) <- label_df$gene_name

  # ============================================================================
  # plot
  # ============================================================================
  if(is.null(sample_group_info)){
    facet_var = sample~gene_name
  }else{
    facet_var = group + sample~gene_name
  }

  # draw track
  if(plot_type == "translatome"){
    p <- ggplot() +
      # geom_col(data = signal_data,
      #          mapping = aes(x = transpos,y = density,
      #                        fill = type),
      #          width = 1) +
      do.call(geom_col,modifyList(list(data = signal_data,
                                       mapping = aes(x = transpos,y = density,
                                                     fill = type),
                                       width = 1),
                                  geom_col_params)) +
      scale_fill_manual(values = signal_col,name = "") +
      xlab("Ribosome Along Transcript (nt)") +
      ylab("Footprint Occupancy")
  }else if(plot_type == "interactome"){
    p <- ggplot() +
      # geom_line(data = signal_data,
      #           mapping = aes(x = codon_pos,y = density,color = gene_name)) +
      do.call(geom_line,modifyList(list(data = signal_data,
                                        mapping = aes(x = codon_pos,y = density,color = gene_name)),
                                   geom_line_params)) +
      # geom_hline(yintercept = 1.5,lty = 'dashed',color = 'red',size = 0.75) +
      do.call(geom_hline,modifyList(list(yintercept = 1.5,lty = 'dashed',
                                         color = 'red',size = 0.75),
                                    geom_hline_params)) +
      # geom_ribbon(data = signal_data,
      #             mapping = aes(x = codon_pos,y = density,
      #                           ymin = density - density_sd,
      #                           ymax = density + density_sd,
      #                           fill = gene_name),
      #             alpha = 0.4) +
      do.call(geom_ribbon,modifyList(list(data = signal_data,
                                          mapping = aes(x = codon_pos,y = density,
                                                        ymin = density - density_sd,
                                                        ymax = density + density_sd,
                                                        fill = gene_name),
                                          alpha = 0.4),
                                     geom_ribbon_params)) +
      line_col_layer +
      ylab('Mean enrichment [AU] (co-IP/total)') +
      xlab('Ribosome position \n (codons/amino acids)')
  }

  # whether add gene structure
  if(add_gene_struc == TRUE){
    ptheme <- p +
      ggnewscale::new_scale_fill() +
      geom_rect(data = structure_df,
                aes(xmin = start,xmax = end,ymin = ymin,ymax = ymax,
                    fill = region),show.legend = FALSE) +
      scale_fill_manual(values = structure_col,name = "") +
      geom_text(data = structure_df,
                mapping = aes(x = (start + end)/2,y = structure_label_y,label = region),
                size = structure_label_size,fontface = "bold.italic")
  }else{
    ptheme <- p
  }

  # add theme details
  ptheme <- ptheme +
    theme_bw() +
    ggh4x::facet_nested(facet_var,
                        nest_line = element_line(linetype = "solid"),
                        scales = "free",independent = "y",
                        switch = "y",
                        labeller = labeller(gene_name = new_label)) +
    ggh4x::facetted_pos_scales(y = range_layer) +
    jj_theme() +
    ggh4x::force_panelsizes(rows = rep(c(rep(panel_size[1],n_sample),panel_size[2]),
                                       length(unique(signal_data$gene_name))))

  # whether show y ticks
  if(show_y_ticks == TRUE){
    ptheme <- ptheme +
      theme(strip.background = element_rect(fill = NA,color = NA),
            strip.text.y.left = element_text(angle = 0),
            # axis.text.x = element_text(),
            # axis.ticks.x = element_line(),
            axis.text.y = element_text(),
            axis.ticks.y = element_line(),
            legend.text = element_text(face = "bold"),
            panel.background = element_rect(fill = alpha(background_col,0.5))) +
      coord_cartesian(clip = "off")
  }else{
    ptheme <- ptheme +
      zplyr::geom_abs_text(data = track_range,
                           aes(xpos = range_pos[1],ypos = range_pos[2],
                               label = rg_label)) +
      theme(strip.background = element_rect(fill = NA,color = NA),
            strip.text.y.left = element_text(angle = 0),
            # axis.text.x = element_text(),
            # axis.ticks.x = element_line(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.text = element_text(face = "bold"),
            panel.background = element_rect(fill = alpha(background_col,0.5))) +
      coord_cartesian(clip = "off")
  }

  # whether show x ticks
  if(show_x_ticks == TRUE){
    ptheme <- ptheme +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line())
  }else{
    ptheme <- ptheme +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }

  # ============================================================================
  # remove panel borders
  # ============================================================================
  if(remove_all_panel_border == TRUE | remove_trans_panel_border == TRUE){
    g <- ggplotGrob(ptheme)

    # remove_all_panel_border = T
    col_num <- length(unique(signal_data$gene_name))

    # defnie panels
    # remove_trans_panel_border = T
    if(remove_all_panel_border == TRUE){
      panel_num <- 1:(col_num*(n_sample + 1))
    }else if(remove_trans_panel_border == TRUE){
      start = col_num*n_sample + 1
      end = start + col_num - 1
      panel_num <- start:end
    }else{
      message("Should not be both TRUE!")
    }

    # grobs editting
    for (i in panel_num + 1) {
      grobs_border <- grid::grid.ls(g$grobs[[i]],print = FALSE)
      panel_boder_name <- grobs_border[["name"]][length(grobs_border[["name"]])]
      g$grobs[[i]]$children[[panel_boder_name]]$gp$col <- NA
    }
    grid::grid.newpage()
    grid::grid.draw(g)
  }else{
    return(ptheme)
  }

}



#' Plot mapping information from mapinfo files
#'
#' This function reads in mapinfo files and extracts mapping information for each file.
#' The extracted information is then plotted either as a barplot or a table.
#'
#' @param mapinfo_file A character vector of file paths to the mapinfo files.
#' @param file_name An optional character vector of names to assign to the samples.
#' If NULL, the sample name will be extracted from the mapinfo file paths.
#' @param plot_type A string indicating the type of plot to generate.
#' Either "barplot" for a stacked barplot of mapping percentages, or "table" for
#' a table of mapping counts.
#' @param geom_col_params A list of parameters to pass to the \code{geom_col} function.
#' Useful for modifying the appearance of the barplot.
#'
#' @return A barplot or table depending on the \code{plot_type} parameter.
#'
#' @examples
#' \dontrun{
#' # Generate a barplot of mapping percentages
#' plot_mapinfo(mapinfo_file = c("sample1_mapinfo.txt", "sample2_mapinfo.txt"), plot_type = "barplot")
#'
#' # Generate a table of mapping counts
#' plot_mapinfo(mapinfo_file = c("sample1_mapinfo.txt", "sample2_mapinfo.txt"), plot_type = "table")
#' }
#'
#' @export
plot_mapinfo <- function(mapinfo_file = NULL,file_name = NULL,
                         plot_type = c("barplot","table"),
                         geom_col_params = list()){
  plot_type <- match.arg(plot_type,c("barplot","table"))

  # define extract func
  extrac_fun <- function(x){
    as.numeric(sapply(strsplit(x,split = "\\("), "[",1))
  }

  # assign file name
  if(is.null(file_name)){
    sample_name <- sapply(strsplit(mapinfo_file,split = "mapinfo.txt"), "[",1)
  }else{
    sample_name <- file_name
  }

  # loop for read and extract
  # x = 1
  plyr::ldply(seq_along(mapinfo_file),function(x){
    tmp <- read.delim(mapinfo_file[x],check.names = FALSE)

    lapply(1:4, function(x){
      extrac_fun(tmp[x,])
    }) |> unlist() -> map_reads

    res <- data.frame(sample = sample_name[x],
                      total_mapped = map_reads[1],
                      un_mapped = map_reads[2],
                      uniq_mapped = map_reads[3],
                      multi_mapped = map_reads[4])

    return(res)
  }) -> all_map_df

  # wide to long
  df_long <- reshape2::melt(all_map_df[,-2],id.vars = "sample",
                            variable.name = "map_type",value.name = "reads")

  # plot
  barplot <-
    ggplot(df_long) +
    do.call(geom_col,modifyList(
      list(mapping = aes(x = reads,y = sample,fill = map_type),
           position = position_fill()),geom_col_params)) +
    scale_fill_brewer(palette = "Paired") +
    scale_x_continuous(labels = scales::label_percent()) +
    theme_bw() + xlab("reads percent") +
    jj_theme() + ylab("")

  # return
  if(plot_type == "barplot"){
    barplot
  }else{
    tab <- gridExtra::tableGrob(all_map_df)
    grid::grid.newpage()
    grid::grid.draw(tab)
  }
}
