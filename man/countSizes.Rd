\name{countSizes}
\alias{countSizes}
\title{
Calculate Region Sizes for All Transcripts
}
\description{
Count sizes for whole transcript,, 5'UTR, cds, 3'UTR regions of each transcript.
}
\usage{
countSizes(insets, trxBed)
}
\arguments{
  \item{insets}{
  List of integer, insets in nucleotides for avoiding star tand stop codons.
}
  \item{trxBed}{
  GRanges list for all transcripts.
}
}
\value{
List of sizes of whole transcript, cds, 5-UTR, and 3-UTR for each transcript.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
