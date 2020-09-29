\name{getBamFlagStat}
\alias{getBamFlagStat}
\title{
Extract Flags from A BAM File
}
\description{
Scan bam flags to get statistics of the reads. This requires samtools available from system (either path to samtools is included in user's PATH variable or the module has been loaded in HPC system). 
}
\usage{
getBamFlagStat(bamFiles, outFile)
}
\arguments{
  \item{bamFiles}{
  Character vector, names of bam files.
}
  \item{outFile}{
  Character vector, names of output file.
}
}
\value{
None. Write the flags to output file.
}
\references{
http://www.htslib.org/doc/samtools.html
}
\keyword{ methods }
