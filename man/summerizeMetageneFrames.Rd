\name{summerizeMetageneFrames}
\alias{summerizeMetageneFrames}
\title{
Summerize Metagene Framse from Metagene Position Profile
}
\description{
Summrize reading frames for each metagene position (position of read start relative to cds start or to cds stop position)
}
\usage{
summerizeMetageneFrames(posProfile, pos_col=1)
}
\arguments{
  \item{posProfile}{
  A data frame with three columns for to_start, to_end, and read_len for each reads.
}
  \item{pos_col}{
  Positive integer, column number in the data frame above.
}
}
\value{
A data frame of 3 columns for total counts, frame, metagene position, and rows for each metagene position
}
\keyword{ methods }
