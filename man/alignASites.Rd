\name{alignASites}
\alias{alignASites}
\title{
Align a-site for GAlignments Object
}
\description{
Covert GAligments objects for a transcript to a GAlignments object having the A site nucleotides only. 
}
\usage{
alignASites(asiteOffsets, alns)
}
\arguments{
  \item{asiteOffsets}{
  A data frame of one column, row names are read lengths and column is A sites for each read length.
}
  \item{alns}{
  A GAlignments objects for a transript.
}
}
\value{
A GAlignments object for A site nucleotides only.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }

