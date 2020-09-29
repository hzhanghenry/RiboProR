\name{getAnnotationFromBedFile}
\alias{getAnnotationFromBedFile}
\title{
Read Annotation Information from BED File
}
\description{
Import annotation information from BED file and convert them to two GRanges lists.
}
\usage{
getAnnotationFromBedFile(bed_file)
}
\arguments{
  \item{bed_file}{
Character vector, file name (and path) for annotations in BED format. The file must be in Bed12 format (with thick columns defined).
}
}
\details{
import() function from rtracklayer package is used to bring in annotation informaton to GRanges list.
}
\value{
List of two GRanges objects.  The first one contains all information in the bed file and the second one has cds information only.
}
\references{
https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
}
\keyword{ methods }
