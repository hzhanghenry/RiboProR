\name{getAlignmentsFromBamFile}
\alias{getAlignmentsFromBamFile}
\title{
Read All Reads from BAM File
}
\description{
Scan a BAM file and hold all reads with GRanges list.
}
\usage{
getAlignmentsFromBamFile(bam_file)
}
\arguments{
  \item{bam_file}{
   Character vector, a bam file name (and path).
}
}
\value{
Granges list containing all reads in bam file. Keep the reads unfiltered but hold seqname, start and end, strand information only.
}
\references{
http://samtools.github.io/hts-specs/SAMv1.pdf
}
\keyword{ methods }


