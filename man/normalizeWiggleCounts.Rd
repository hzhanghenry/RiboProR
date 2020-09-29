\name{normalizeWiggleCounts}
\alias{normalizeWiggleCounts}
\title{
Normalize Wiggle Track Counts
}
\description{
Normalize wiggle track counts to a defined total counts.
}
\usage{
normalizeWiggleCounts(all_counts, normalize_factor = NULL)
}
\arguments{
  \item{all_counts}{
  List of list, each list element has two numeric vectors for read counts forward and reverse strand at base pair level.
}
  \item{normalize_factor}{
  positive numeric, scaling factor. If not provided, total wiggle counts will be scaled to 1000000000.
}
}
\value{
List of list, normalized counts for each strand.
}
\references{
https://github.com/ingolia-lab/RiboSeq
}
\keyword{ methods }
