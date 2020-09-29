\name{transcriptGRanges}
\alias{transcriptGRanges}
\title{
Get GRanges from BED Content for Exons of A Transcript
}
\description{
Extract exon information from a transcript GRanges object and hold with GRanges list.
}
\usage{
transcriptGRanges(bedGRange)
}
\arguments{
  \item{bedGRange}{
  A GRange object for a transcript with blocks in metadata columns.
}
}
\value{
GRanges list with one GRange per exon.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
