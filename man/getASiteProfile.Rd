\name{getASiteProfile}
\alias{getASiteProfile}
\title{
Find All A-site for Alignments of A Transcript
}
\description{
Count number of aligned a-sites at each base position of a transcripts.
}
\usage{
getASiteProfile(asiteOffsets, bamfile, trx)
}
\arguments{
  \item{asiteOffsets}{
  A data frame of one column with read length as row names and having a-sites for each read length in the column.
}
  \item{bamfile}{
  Character vector, name of a bam file (and path).
}
  \item{trx}{
  GRange object for a transcript.
}
}
\value{
Integer vector,  total aligned a-site at each base position of a transcript.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }

