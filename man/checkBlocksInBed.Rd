\name{checkBlocksInBed}
\alias{checkBlocksInBed}
\title{
Validate the Information of Blocks Defined in BED File
}
\description{
Check out if the block information defined in bed file are correct.
Block count must be positvie integer. Block sizes and block starts must be comma seperated list.
Number of block sized and starts must match to the block counts.
}
\usage{
checkBlocksInBed(bed_info)
}
\arguments{
  \item{bed_info}{
A data frame with 12 columns read from bigBed file.It should be checked before call this function and must have required standard column names.
}
}
\value{
Logic, TRUE for everything is correct. Otherwise FALSE.
}
\references{
https://genome.ucsc.edu/FAQ/FAQformat#format1.7
}
\keyword{ methods }
