globalVariables(c("counts", "dfqc", "feature", "framest", "group", "num",
                  "normaledReads",":=","density", "gene_name", "if_else",
                  "max_v", "min_v", "remove_chrom_panel_border", "rg_label",
                  "trans_id", "transpos", "type", "ymax", "ymin",".","region"))

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
  # theme setting
  jj_theme <-
    theme(strip.background = element_rect(),
          axis.text = element_text(color = "black",size = rel(1)),
          axis.title = element_text(color = "black",size = rel(1.2)),
          strip.text = element_text(color = "black",size = rel(1),face = "bold"),
          panel.grid = element_blank())

  # ===================================================================
  # plot layers
  # ===================================================================
  # plot
  pmain <-
    ggplot() +
    theme_bw() +
    jj_theme +
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
                              fill = "#8294C4",width = 0.75),
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
#' @param shift A numeric value specifying the p site to shift, default 12.
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
                              shift = 12,
                              geom_col_list = list(),
                              facet_wrap_list = list()){
  options(warn=-1)
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)
  # ===================================================================
  # theme setting
  jj_theme <-
    theme(strip.background = element_rect(),
          axis.text = element_text(color = "black",size = rel(1)),
          axis.title = element_text(color = "black",size = rel(1.2)),
          strip.text = element_text(color = "black",size = rel(1),face = "bold"),
          panel.grid = element_blank())

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
    jj_theme +
    xlab('Read length') + ylab('Reads numbers(k)') +
    do.call(facet_wrap,
            modifyList(list(~sample,scales = 'free',ncol = 2),
                       facet_wrap_list)) +
    do.call(geom_col,modifyList(
      list(data = df,
           mapping = aes(x = !!vars_f[[1]] + shift,y = normaledReads/1000,
                         fill = factor(!!vars_f[[2]])),
           width = 1,
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
#' @param structure_col Color of the gene structure rectangles.
#' @param background_col Color of the plot background. Default is "grey90".
#' @param range_pos Position of the range labels. Default is c(0.85,0.85).
#' @param signal_col A named vector of colors for ribosome and RNA signals.
#' Default is NULL, which means using a default color scheme.
#' @param panel_size A numeric vector of length 2 specifying the width and height
#' of each panel. Default is NULL, which means using a default panel size.
#' @param remove_all_panel_border Logical indicating whether to remove all panel
#' borders. Default is FALSE.
#' @param remove_trans_panel_border Logical indicating whether to remove the
#' border of the transcript panel only. Default is FALSE.
#' @param gene_order An optional character vector specifying the order of genes
#' in the plot. By default, genes are ordered by their appearance in the input data.
#' @param sample_order An optional character vector specifying the order of
#' samples in the plot. By default, samples are ordered by their appearance in
#' the input data.
#' @param fixed_col_range Logical indicating whether to use a fixed color range
#' for all panels. Default is TRUE.
#' @param show_ribo_only Wthether show only ribo density track. Default is FALSE.
#' @param sample_group_info The groups for samples, giving a named list with
#' samples, default NULL.
#'
#' @import zplyr
#' @import ggh4x
#' @importFrom grid grid.ls grid.draw grid.newpage
#'
#' @export
track_plot <- function(signal_data = NULL,
                       gene_anno = NULL,
                       structure_col = NULL,
                       background_col = "white",
                       range_pos = c(0.85,0.85),
                       signal_col = NULL,
                       panel_size = NULL,
                       remove_all_panel_border = FALSE,
                       remove_trans_panel_border = FALSE,
                       gene_order = NULL,
                       sample_order = NULL,
                       fixed_col_range = TRUE,
                       show_ribo_only = FALSE,
                       sample_group_info = NULL){
  # ==============================================================================
  # process data
  # ==============================================================================
  # show_ribo_only = T
  if(show_ribo_only == TRUE){
    signal_data <- signal_data %>% dplyr::filter(type != "rna")
  }else{
    signal_data <- signal_data %>%
      dplyr::mutate(density = dplyr::if_else(type %in% "rna",-density,density))

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
  # gene strctures
  # ==============================================================================
  # load geneinfo
  geneInfo <- read.table('longest_info.txt')
  colnames(geneInfo) <- c('id','gene_name','gene_id','trans_id','chr','strand',
                          'cds_region','exon_region','utr5','cds','utr3')
  geneInfo <- geneInfo %>%
    dplyr::filter(gene_name %in% unique(signal_data$gene_name))

  # x = 1
  structure_df <- plyr::ldply(1:nrow(geneInfo), function(x){
    tmp <- geneInfo[x,]
    df <- data.frame(gene_name = tmp$gene_name,
                     start = c(0,tmp$utr5,tmp$utr5 + tmp$cds),
                     end = c(tmp$utr5,tmp$utr5 + tmp$cds,tmp$utr5 + tmp$cds + tmp$utr3),
                     ymin = c(-0.5,-1,-0.5),
                     ymax = c(0.5,1,0.5),
                     sample = "trans",
                     region = c("5UTR","CDS","3UTR"),
                     group = NA)

    return(df)
  })

  # ==============================================================================
  # panel range settings
  # ==============================================================================
  n_sample = length(unique(unique(signal_data$sample)))

  # fixed_col_range = F
  if(fixed_col_range == TRUE){
    track_range <- signal_data %>%
      dplyr::group_by(gene_name) %>%
      dplyr::summarise(min_v = min(density),
                       max_v = max(density)) %>%
      dplyr::mutate(rg_label = paste("[",signif(min_v,digits = 2),"-",
                                     signif(max_v,digits = 2),"]",sep = "")) %>%
      replicate(n_sample,.,simplify = F) %>%
      do.call("rbind",.)

    track_range$sample <- rep(unique(signal_data$sample),each = n_sample)
  }else{
    track_range <- signal_data %>%
      dplyr::group_by(gene_name,sample) %>%
      dplyr::summarise(min_v = min(density),
                       max_v = max(density)) %>%
      dplyr::mutate(rg_label = paste("[",signif(min_v,digits = 2),"-",
                                     signif(max_v,digits = 2),"]",sep = ""))
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

  # y axis range layer
  range_layer <- lapply(1:(nrow(track_range) + length(unique(track_range$gene_name))), function(x){
    if(x <= nrow(track_range)){
      tmp <- track_range[x,]
      scale_y_continuous(limits = c(tmp$min_v,tmp$max_v))
    }else{
      scale_y_continuous(limits = c(-2,1))
    }
  })

  # ==============================================================================
  # paras settings
  # ==============================================================================
  # theme setting
  jj_theme <-
    theme(strip.background = element_rect(),
          axis.text = element_text(color = "black",size = rel(1)),
          axis.title = element_text(color = "black",size = rel(1.2),face = "bold"),
          strip.text = element_text(color = "black",size = rel(1),face = "bold"),
          strip.text.x = element_text(color = "black",size = rel(1),face = "bold.italic"),
          panel.grid = element_blank())

  # color
  if(is.null(signal_col)){
    signal_col <- c("ribo" = "#D21312","rna" = "#009FBD")
  }else{
    signal_col <- signal_col
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
  label_df <- signal_data %>%
    dplyr::select(gene_name,trans_id) %>% unique()

  new_label <- paste(label_df$gene_name,label_df$trans_id,sep = "\n")
  names(new_label) <- label_df$gene_name

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
  # ============================================================================
  # plot
  # ============================================================================
  if(is.null(sample_group_info)){
    facet_var = sample~gene_name
  }else{
    facet_var = group + sample~gene_name
  }

  # draw
  p <-
    ggplot() +
    geom_col(data = signal_data,
             mapping = aes(x = transpos,y = density,
                           fill = type),
             width = 1) +
    scale_fill_manual(values = signal_col,name = "") +
    # geom_hline(yintercept = 0,lty = "solid",color = "black") +
    ggnewscale::new_scale_fill() +
    geom_rect(data = structure_df,
              aes(xmin = start,xmax = end,ymin = ymin,ymax = ymax,
                  fill = region),show.legend = FALSE) +
    scale_fill_manual(values = structure_col,name = "") +
    geom_text(data = structure_df,
              mapping = aes(x = (start + end)/2,y = -1.5,label = region),
              size = 3,fontface = "bold.italic") +
    zplyr::geom_abs_text(data = track_range,
                         aes(xpos = range_pos[1],ypos = range_pos[2],
                             label = rg_label)) +
    theme_bw() +
    xlab("Ribosome Along Transcript (nt)") +
    ylab("Footprint Occupancy(upper)\n RNA Covergae(bottom)") +
    ggh4x::facet_nested(facet_var,
                        nest_line = element_line(linetype = "solid"),
                        scales = "free",independent = "y",
                        switch = "y",
                        labeller = labeller(gene_name = new_label)) +
    ggh4x::facetted_pos_scales(y = range_layer) +
    jj_theme +
    theme(strip.background = element_rect(fill = NA,color = NA),
          strip.text.y.left = element_text(angle = 0),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(face = "bold"),
          panel.background = element_rect(fill = alpha(background_col,0.5))) +
    ggh4x::force_panelsizes(rows = rep(c(rep(panel_size[1],n_sample),panel_size[2]),
                                       length(unique(signal_data$gene_name))))

  # ============================================================================
  # remove panel borders
  # ============================================================================
  if(remove_all_panel_border == TRUE | remove_trans_panel_border == TRUE){
    g <- ggplotGrob(p)

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
    return(p)
  }

}
