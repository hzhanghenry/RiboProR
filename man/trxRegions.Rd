\name{trxRegions}
\alias{trxRegions}
\title{
Calculate Transcript Regions for UTR and CDS
}
\description{
Calculate length of whole transcript, 5'UTR, cds, and 3'UTR regions for a transcript.
}
\usage{
trxRegions(trx, cds, insetUtr5Start=0, insetUtr5End=0, insetCdsStart=0, 
	insetCdsEnd=0, insetUtr3Start=0, insetUtr3End=0)
}
\arguments{
  \item{trx}{
  GRanges object for a transcript.
}
  \item{cds}{
  GRanges object for cds in a transcript.
}
  \item{insetUtr5Start}{
  Positive integer, adjustment after start position of 5'UTR.
}
  \item{insetUtr5End}{
  Positive integer, adjustment before end position of 5'UTR .
}
  \item{insetCdsStart}{
 Positive integer, adjustment after start position of cds.
}
  \item{insetCdsEnd}{
 Positive integer, adjustment before end position of cds.
}
  \item{insetUtr3Start}{
 Positive integer, adjustment after start position of 3'UTR.
}
  \item{insetUtr3End}{
  Positive integer, adjustment before end position of 3'UTR .
}
}
\value{
GRanges list representing cds, 5-UTR, and 3-UTR for a transcript.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
