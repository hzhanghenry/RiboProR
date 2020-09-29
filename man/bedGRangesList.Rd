\name{bedGRangesList}
\alias{bedGRangesList}
\title{
Get GRanges List for All Transcripts in A BED File
}
\description{
For each transcript in a BED file GRanges, get the GRanges of its exons as per transcriptGRanges and collec these into a GRangesList with names taken from the name metadata column of the BED file GRanges.
}
\usage{
bedGRangesList(bed)
}
\arguments{
  \item{bed}{
  GRanges list with block information defined.
}
}
\value{
A GRanges list object for exons of all transcripts.
}
\references{
Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379
}
\keyword{ methods }

