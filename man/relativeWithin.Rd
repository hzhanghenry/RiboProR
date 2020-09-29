\name{relativeWithin}
\alias{relativeWithin}
\title{
Convert a Genomic Position to Position Relative to CDS Start
}
\description{
Convert an absolute genomic query position (qpos) to a transcript-relative coordinate position within a transcript.
}
\usage{
relativeWithin(qpos, outer)
}
\arguments{
  \item{qpos}{
  A positive integer, a genomic coordinate inside of a transcript.
}
  \item{outer}{
  GRanges object, genomic positions of a transctipt.
}
}
\value{
Positive integer, a poistion relative to the start of transctipt or NA if the absolute position is not in any exon (no hit).
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
