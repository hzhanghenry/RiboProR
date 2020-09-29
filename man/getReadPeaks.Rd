\name{getReadPeaks}
\alias{getReadPeaks}
\title{
Find Read Peaks from Metagene Table
}
\description{
Find count peaks for each read length from metagene length table.
}
\usage{
getReadPeaks(metagene, min_range, max_range)
}
\arguments{
  \item{metagene}{
  A data frame with columns for each read length and rows for counts in each position of metagene.
}
  \item{min_range}{
  Positive integer, minimum index of the row to find peak.
}
  \item{max_range}{
  Positive integer, maximum index of the row to find peak.
}
}
\value{
An integet vector, distance from metagene start for each read length.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
