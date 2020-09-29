\name{getWiggleCounts}
\alias{getWiggleCounts}
\title{
Get A-site Counts for Wiggle File Generation
}
\description{
Count total number of aligned a-sites for full genome or transcripts only.
}
\usage{
getWiggleCounts(bam_file, asite_table, annot_bed_file)
}
\arguments{
  \item{bam_file}{
  Chraracter vector, name of bam file (and path) 
}
  \item{asite_table}{
  A data frame with 1 column for a-site and rownames for read length. A file name is also accepted.
}
  \item{annot_bed_file}{
  Chraracter vector, name of a bed file for gene/transcript annotation.
}
}
\value{
List of list with each sub-list is wiggle counts for forward and reverse strand of one chromosome.
}
\references{
https://github.com/ingolia-lab/RiboSeq
}
\note{
This function counts coverage for forward and reverse strand separately and the coverage is on the A-site only.
}
\keyword{ methods }


