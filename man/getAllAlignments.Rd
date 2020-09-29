\name{getAllAlignments}
\alias{getAllAlignments}
\title{
Get All Alignments for A Transcript from BAM File
}
\description{
Get all reads overlapping the genomic extent of a primary transcript on the correct strand, including unspliced and purely intronic reads.
}
\usage{
getAllAlignments(bamfile, trx)
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

