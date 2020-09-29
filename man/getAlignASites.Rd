\name{getAlignASites}
\alias{getAlignASites}
\title{
Get Aligned A sites for A Transcript
}
\description{
Scan  bam file for all reads of a transcript then get aligned a-sites for the transcript.
}
\usage{
getAlignASites(asiteOffsets, bamfile, trx)
}
\arguments{
  \item{asiteOffsets}{
  A data frame of one column with read length as row names and a-sites for each read length in column.
}
  \item{bamfile}{
  Character vector, name of a bam file (and path).
}
  \item{trx}{
  GRange object for a transcript.
}
}
\value{
A GAlignments object for A site nucleotides only.
}
\references{
Oroginal code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
