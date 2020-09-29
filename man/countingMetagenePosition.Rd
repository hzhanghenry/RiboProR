\name{countingMetagenePosition}
\alias{countingMetagenePosition}
\title{
Count Metagene Positions for Aligments
}
\description{
Find the positions relative to cds start and stop positions for all alignments.
}
\usage{
countingMetagenePosition(alignments, annotation, shift_start=0, shift_end=0)
}
\arguments{
  \item{alignments}{
  GRanges list for short reads in a bam file.
}
  \item{annotation}{
  GRanges list of cds annotations.
}
  \item{shift_start}{
  Positive integer, shift this distance to 5'end from start position of cds.
}
  \item{shift_end}{
  Positive integer, shift this distance to 3'end from end position of cds.
}
}
\value{
A data frame with 3 columns for 3 columns: distance to cds start, distance to cds end, and read length, for each aligments.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
