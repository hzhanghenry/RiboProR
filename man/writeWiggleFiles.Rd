\name{writeWiggleFiles}
\alias{writeWiggleFiles}
\title{
Write Aligned A-site Counts to Wiggle Files
}
\description{
Write aligned a-site Counts on each chromosome to wiggle track files by strand.
}
\usage{
writeWiggleFiles(all_counts, bam_file, is.normalized=FALSE)
}
\arguments{
  \item{all_counts}{
  List of list, each list element has two numeric vectors for read counts on forward and reverse strand at base pair level.
}
  \item{bam_file}{
  character vector, name of bam file from which the wiggle track data is generated.
}
  \item{is.normalized}{
  Logic, if the all_counts are normalized.
}
}
\value{
None. Write files only
}
\references{
https://github.com/ingolia-lab/RiboSeq
https://genome.ucsc.edu/goldenpath/help/wiggle.html
}
\keyword{ methods }

