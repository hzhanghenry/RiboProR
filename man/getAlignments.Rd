\name{getAlignments}
\alias{getAlignments}
\title{
Get Alignments for A TRanscript from BAM File
}
\description{
Read in all reads compatible with a processed transcript on the correct strand, as per findSpliceOverlaps.
}
\usage{
getAlignments(bamfile, trx)
}
\arguments{
  \item{bamfile}{
  Character vector, the name (and path) of the bam file to be read.
}
  \item{trx}{
  GRanges object for a transcript.
}
}
\value{
A GAlignments object with all reads overlapping with the transcript.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }

