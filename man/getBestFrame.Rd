\name{getBestFrame}
\alias{getBestFrame}
\title{
Find the Best Reading Frame
}
\description{
Find the best frame from Frame 0~2 for each read length.
}
\usage{
getBestFrame(frame_table)
}
\arguments{
  \item{frame_table}{
  An data frame with counts of 3 readind frames for each read length.
}
}
\value{
Positive integer vector, frame number for each read length,
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
