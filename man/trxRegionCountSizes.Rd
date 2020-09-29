\name{trxRegionCountSizes}
\alias{trxRegionCountSizes}
\title{
Get Region Sized of A Transcript
}
\description{
Calculate sizes of whole transcript, cds, 5'UTR, and 3'UTR of a transcript defined in bed file.
}
\usage{
trxRegionCountSizes(insets, trx, cds)
}
\arguments{
  \item{insets}{
  List of integers of length 2, insets in nucleotides to avoid start and stop positions. 
}
  \item{trx}{
  GRange object for transcripts.
}
  \item{cds}{
  IRanges object for cds of the transcripts.
}
}
\value{
List of integers pf length 4 for sizes of transcript, cds, 5-UTR, and 3-UTR..
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
