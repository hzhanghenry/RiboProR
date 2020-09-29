\name{initializeMetageneTable}
\alias{initializeMetageneTable}
\title{
Initialize A New Metagene Table
}
\description{
Initialize an empty matrix with rows for base positions of metagene and columns for read length.
}
\usage{
initializeMetageneTable(meta_start=-100, meta_end=100, min_len=25, max_len=34)
}
\arguments{
  \item{meta_start}{
  Negative integer, distance before cds start for metagene start.
}
  \item{meta_end}{
  Positive integer, distance after cds end for metagene end.
}
  \item{min_len}{
  Positive integer, minimum length of a read.
}
  \item{max_len}{
  Positive integer, maximum length of a read.
}
}
\value{
A matrix of 0s(zero) with rows for metagene positions and columns for read length.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
