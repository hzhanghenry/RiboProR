\name{getCDSGRanges}
\alias{getCDSGRanges}
\title{
Extract CDS Regions from Transcript Definition
}
\description{
Covert GRanges of transcripts to GRanges of CDS regions.
}
\usage{
getCDSGRanges(bed_file)
}
\arguments{
  \item{bed_file}{
character vector, name (and path) of a bed file defining transcripts and relevant CDS regions.
}
}
\value{
GRange list for CDS regions.
}
\keyword{ methods }
