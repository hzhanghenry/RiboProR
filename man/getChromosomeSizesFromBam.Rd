\name{getChromosomeSizesFromBam}
\alias{getChromosomeSizesFromBam}
\title{
Get Chromosome Sizes from A BAM File
}
\description{
Extract length of each chromosome from bam file header.
}
\usage{
getChromosomeSizesFromBam(bam_file)
}
\arguments{
  \item{bam_file}{
  Character vector, bam file name (and path).
}
}
\value{
A data frame with 2 columns for chromosome names and lengths.
}
\keyword{ methods }
