% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilsFuncs.R
\name{fetch_trRNA_from_NCBI}
\alias{fetch_trRNA_from_NCBI}
\title{Fetch tRNA and rRNA sequences from NCBI database}
\usage{
fetch_trRNA_from_NCBI(species = "Homo sapiens", output_file = NULL)
}
\arguments{
\item{species}{A character string specifying the scientific name of the target
species. Defaults to "Homo sapiens".}

\item{output_file}{A character string specifying the name of the output file.
If not provided, the default name is "<species>_trRNA.fa".}
}
\value{
A FASTA-formatted file containing the tRNA and rRNA sequences for the
specified species is created.
}
\description{
This function allows users to fetch tRNA and rRNA sequences for a given
species from the NCBI nucleotide database.
}
\examples{
\dontrun{
# Fetch tRNA and rRNA sequences for human
fetch_trRNA_from_NCBI("Homo sapiens")

# Fetch tRNA and rRNA sequences for mouse and save to a custom file name
fetch_trRNA_from_NCBI("Mus musculus", "my_mouse_trRNA.fa")
}

}
\references{
\url{https://www.ncbi.nlm.nih.gov/}
}
