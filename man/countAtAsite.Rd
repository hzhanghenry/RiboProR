\name{countAtAsite}
\alias{countAtAsite}
\title{
Count A-site Along One Chromosome
}
\description{
Count number of reads which overlap at A-site on transcripts of one chromosome.
}
\usage{
countAtAsite(alignments, chrom_len, asite_table)
}
\arguments{
  \item{alignments}{
  GRange list for alignments on one chromosome.
}
  \item{chrom_len}{
  Positive integer, length of the chromosome.
}
  \item{asite_table}{
  A data frame with one column for A-sites and rownames for qualified read length.
}

}
\value{
A list of 2 integer vectors for counts on forward and reverse strand at each base pair position.
}
\references{
https://github.com/ingolia-lab/RiboSeq
}
\keyword{ methods }

