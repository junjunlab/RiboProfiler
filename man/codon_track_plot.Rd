% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codon_track_plot.R
\name{codon_track_plot}
\alias{codon_track_plot}
\alias{codon_track_plot,ribosomeObj-method}
\title{Generate Codon Track Plot}
\usage{
codon_track_plot(
  object,
  codon_exp_file = NULL,
  sample_name = NULL,
  group_name = NULL,
  select_gene = NULL,
  cds_region_anno = NULL,
  cds_col = "grey",
  cds_anno_col = "#996600",
  cds_region_width = 5,
  cds_region_scale_x = NULL,
  geom_col_list = list(),
  ...
)

\S4method{codon_track_plot}{ribosomeObj}(
  object,
  codon_exp_file = NULL,
  sample_name = NULL,
  group_name = NULL,
  select_gene = NULL,
  cds_region_anno = NULL,
  cds_col = "grey",
  cds_anno_col = "#996600",
  cds_region_width = 5,
  cds_region_scale_x = NULL,
  geom_col_list = list(),
  ...
)
}
\arguments{
\item{object}{ribosomeObj object.}

\item{codon_exp_file}{A vector of file paths to codon expression data files.}

\item{sample_name}{A vector of sample names corresponding to each codon expression file.}

\item{group_name}{An optional vector of group names for each sample.}

\item{select_gene}{A vector of gene names to be plotted.}

\item{cds_region_anno}{A data frame containing annotation data for the coding sequence (CDS) region,
which must has gene_name,start,end three columns. Default is \code{NULL}.}

\item{cds_col}{A character string specifying the color for the CDS region. Default is \code{"grey"}.}

\item{cds_anno_col}{A character string specifying the color for the CDS annotation. Default is \code{"#996600"}.}

\item{cds_region_width}{A numeric value specifying the width of the CDS region in the plot. Default is \code{5}.}

\item{cds_region_scale_x}{A numeric value specifying the scale of the x-axis for the CDS region. Default is \code{NULL}.}

\item{geom_col_list}{A list of additional parameters to be passed to geom_col().}

\item{...}{Useless args.}
}
\value{
A ggplot object representing the combined track plot for all selected genes.
}
\description{
This function creates a track plot of codon-level ribosome footprint density for selected genes.
}
\details{
This function performs the following steps:
  1. Loads codon expression data from multiple files.
  2. Processes gene annotation data.
  3. Filters and combines data for selected genes.
  4. Creates individual track plots for each gene.
  5. Combines all plots into a single ggplot object.

The resulting plot shows ribosome footprint density along transcript positions for each selected gene,
with separate tracks for each sample or group.
}
\examples{
\dontrun{codon_track_plot(
  object = ribosomeObj,
  codon_exp_file = c("sample1_codon_exp.txt", "sample2_codon_exp.txt"),
  sample_name = c("Sample1", "Sample2"),
  select_gene = c("GENE1", "GENE2"),
  geom_col_list = list(width = 0.8, position = "dodge")
)}

}
