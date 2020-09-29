\name{getFrameTable}
\alias{getFrameTable}
\title{
Generate Reading Frame Table from Frame Profile Table
}
\description{
Convert readind frame profile to readind frame table for total number and fraction of each frame.
}
\usage{
getFrameTable(frameProfile)
}
\arguments{
  \item{frameProfile}{
  A data frame with two columns for frames and length of each read
}
}
\value{
A data frame with 7 columns for read_len, fraction, counts of the three reading frames.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
