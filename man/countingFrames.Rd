\name{countingFrames}
\alias{countingFrames}
\title{
Count Total Number of Reading Frames for All Aligments
}
\description{
Find the reading frames based on the start position of alignments.
}
\usage{
countingFrames(alignments, cds_ranges, inset_5=34, inset_3=31)
}
\arguments{
  \item{alignments}{
  GRanges list for all reads from a bam file.
}
  \item{cds_ranges}{
  GRanges list for all cds annotations.
}
  \item{inset_5}{
  positive integer, shift distance after start position of cds to avoid start codon.
}
  \item{inset_3}{
  positive integer, shift distance before end position of cds to avoid stop codon.
}
}
\value{
A data frame of two columns for reading frames and read length of each read.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
