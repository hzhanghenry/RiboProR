\name{getPositionProfileTable}
\alias{getPositionProfileTable}
\title{
Generate Metagene Table
}
\description{
Convert position profile table to metagene table.
}
\usage{
getPositionProfileTable(positionProfile=NULL, at_which=1, 
	meta_start=-100, meta_end=100)
}
\arguments{
  \item{positionProfile}{
  A data frame of 3 columns for relative position to cds start, to cds end, and read length.
}
  \item{at_which}{
  A positive integer,  1 for to start position or 2 for to end position.
}
  \item{meta_start}{
  Negative integer, position before cds start in metagene table 
}
  \item{meta_end}{
  Positive integer, position after cds end inmetagene table
}
}
\value{
A numeric matrix with rows for metagene positions and columns for read length.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
