% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilsFuncs.R
\name{trRNA_index_build}
\alias{trRNA_index_build}
\title{Build an index for tRNA and rRNA sequences using Bowtie2.}
\usage{
trRNA_index_build(trRNA_file = NULL, prefix = NULL, threads = 8)
}
\arguments{
\item{trRNA_file}{A character string specifying the path to the FASTA file
containing tRNA and rRNA sequences.}

\item{prefix}{A character string specifying the prefix for the output index
files.}

\item{threads}{An integer specifying the number of threads to use for building
the index.}
}
\value{
Nothing is returned; the function saves the index files in the
"0.index-data/rtRNA-index/" directory.
}
\description{
This function builds an index for tRNA and rRNA sequences using Bowtie2.
The resulting index can be used for mapping reads against tRNA/rRNA sequences
to remove them from downstream analysis, as they can interfere with gene
expression measurement.
}
\examples{
\dontrun{
trRNA_index_build("path/to/trRNA.fasta", "trRNA_bt2_index", 8)
}

}
