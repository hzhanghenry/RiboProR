\name{transcriptCdsIRanges}
\alias{transcriptCdsIRanges}
\title{
Convert Exons in GRanges Object to IRanges Object
}
\description{
Convert the thick block information in a GRanges object to transcript-relative IRanges object.
}
\usage{
transcriptCdsIRanges(trx, thickStart, thickEnd)
}
\arguments{
  \item{trx}{
  GRanges object for a transcript with block information.
}
  \item{thickStart}{
  Positive integer, start position of CDS in the transcript.
}
  \item{thickEnd}{
  Positive integer, end position of CDS in the transcript.
}
}
\value{
An IRanges object for relative start and end position of CDS.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
