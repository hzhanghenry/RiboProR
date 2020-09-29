\name{absoluteOuter}
\alias{absoluteOuter}
\title{
Convert A Relative Position to An Absolute Genomic Coordinate
}
\description{
Convert a transcript-relative position to an absolute genomic coordinate based on a transcript whose coordinates are given by a GRanges of exons (outer).
}
\usage{
absoluteOuter(qpos, outer)
}
\arguments{
  \item{qpos}{
  positive integer, a position relative to the start of a transcript.
}
  \item{outer}{
  GRanges object, genomic positions of a transctipt.
}
}
\value{
Positive integer, the geneomic coordinate of the transcript-relative position.
}
\references{
Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379
}
\keyword{ methods }