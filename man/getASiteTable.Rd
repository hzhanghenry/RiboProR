\name{getASiteTable}
\alias{getASiteTable}
\title{
Generate A-site Table
}
\description{
Generate a-site table from the three position profile tables.
}
\usage{
getASiteTable(frame_table, at_start, at_end, parameters)
}
\arguments{
  \item{frame_table}{
  A data frame with 7 columns for read_len, fraction, counts of the three reading frames.
}
  \item{at_start}{
  A numeric matrix with rows for metagene positions realtive to cds start and columns for read length.
}
  \item{at_end}{
  A numeric matrix with rows for metagene positions realtive to cds end and columns for read length.
}
  \item{parameters}{
  List of length 13, all parameters for a-site calculation.
}
}
\value{
A data frame with rows for each read length and columns (total of 9) for read length, fraction of each length, peak relative to cds start, peak relative to cds end, other information, fraction of frame0, frame1, frame2, and best a-site.
}
\references{
(https://github.com/ingolia-lab/RiboSeq)
}
\keyword{ methods }
