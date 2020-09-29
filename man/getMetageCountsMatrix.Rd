\name{getMetageCountsMatrix}
\alias{getMetageCountsMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Matrix of Reads Counts for Metagene Positions
}
\description{
By giving a list of metagene frame table files, generate a matrix for metagene counts plot.
Rows of the matrix are reads counts on each defined metagene position and columns are samples. 
The matrix is sorted by maximum value of each column in decreasing order so that the y range 
of the plot will be automatically decided.
}
\usage{
getMetageCountsMatrix(count_files, from_position, to_position)
}
\arguments{
  \item{count_files}{
  Character vector, names (and path) of metagene reads count files.
}
  \item{from_position}{
  Integer, start position of metagene to plot.
}
  \item{to_position}{
  Integer, stop position of metagene to plot 
}
}
\value{
Data matrix with read counts along metagene position from multiple samples.
}
\author{
Henry Zhang
}
\keyword{ methods }
