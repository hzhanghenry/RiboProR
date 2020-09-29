\name{irangeOuter}
\alias{irangeOuter}
\title{
Get Genomic Ranges for Transcript-Relative IRanges
}
\description{
Convert transcript-relative query IRanges to an absolute genomic GRanges based on a transcript coordinates in a GRanges of exons.
}
\usage{
irangeOuter(qranges, outer)
}
\arguments{
  \item{qranges}{
  IRanges object with relative positions to a start of transcript.
}
  \item{outer}{
  GRanges object, genomic positions of a transctipt.
}
}
\value{
A GRanges object for the genomic positions.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
