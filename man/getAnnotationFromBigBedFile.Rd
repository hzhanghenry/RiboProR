\name{getAnnotationFromBigBedFile}
\alias{getAnnotationFromBigBedFile}
\title{
Extract Annotation Information from BigBed File
}
\description{
Read bigBed file (12 columns) into data frame and check the contents.
}
\usage{
getAnnotationFromBigBedFile(bed_file, has.header=FALSE, sepcial_chrom=NULL)
}
\arguments{
  \item{bed_file}{
Character vector, name (and path) of the bigbed file.
}
  \item{has.header}{
Logic, if the bigbed file has column headers. Default is FALSE.
}
  \item{sepcial_chrom}{
character vector, chromosome names other than digits, roman numbers, and "X", "Y", "M".
}
}
\details{
Fields in bigBed file: 
1)  chrom:      chromosome name of each feature
2)  chromStart: start position of each feature
3)  chromEnd:   end position of each feature
4)  name:       gene name of each feature
5)  score:      used for graphic display only
6)  strand:     chromosome strand of each feature
7)  thickStart: starting position drawn thickly
8)  thickEnd:	 ending position drawn thickly
9)  itemRgb:    used for graphic display only
10) blockCount: number of blocks of each feature
11) blockSizes: comma-separated list of block sizes
12) blokStarts: comma-separated list of block starts, relative to chromStart

chromStart and chromEnd are 0-based and half-open.
}
\value{
A data frame with 12 columns.
}
\references{
https://genome.ucsc.edu/FAQ/FAQformat#format1.7
}
\keyword{ methods }