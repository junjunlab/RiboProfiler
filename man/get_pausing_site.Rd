% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_pausing_site.R
\name{get_pausing_site}
\alias{get_pausing_site}
\alias{get_pausing_site,ribosomeObj-method}
\title{Detect Pausing Sites}
\usage{
get_pausing_site(
  object,
  min_counts = 64,
  norm_type = c("average", "rpm"),
  window = 1000,
  ...
)

\S4method{get_pausing_site}{ribosomeObj}(
  object,
  min_counts = 64,
  norm_type = c("average", "rpm"),
  window = 1000,
  ...
)
}
\arguments{
\item{object}{ribosomeObj object.}

\item{min_counts}{An integer specifying the minimum counts threshold for
considering a site as a pausing site. Default is \code{64}.}

\item{norm_type}{The nomalization methods for ribosome density. "average" is calculated by
the count at each position divided by mean density across cds region. "rpm"
is calculated by the count at each position divided by the total counts and multiplied with 10^6.
Default is "average".}

\item{window}{An integer specifying the window size for detecting pausing
sites. Default is \code{1000}.}

\item{...}{Useless args.}
}
\value{
This function does not return a value. It writes output files to the
"pausing_site_data" directory.
}
\description{
This function processes normalized ribosome profiling data to detect pausing sites.
It writes intermediate normed data files and uses a Python script to perform
the actual detection of pausing sites. This method can be refered to pausepred
[PausePred and Rfeet: webtools for inferring ribosome pauses and visualizing footprint density from ribosome profiling data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6140459/)
}
\examples{
\dontrun{
# Example usage:
get_pausing_site(object = ribosomeObj,
                 min_counts = 64,
                 window = 500)
}
}
