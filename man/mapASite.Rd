\name{mapASite}
\alias{mapASite}
\title{
Mapping A-site for Each Alignment/Read
}
\description{
Convert the GRange object to one base width GRanges based on relevant a-site.
}
\usage{
mapASite(alignments, asite_table)
}
\arguments{
  \item{alignments}{
  GRange object of alignments on one strand to be filtered out.
}
  \item{asite_table}{
  A data frame of 1 columns to hold a-site for  each read length which are represented by row names.
}
}
\value{
GRange object with one base width.
}
\keyword{ methods }

