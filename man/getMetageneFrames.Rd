\name{getMetageneFrames}
\alias{getMetageneFrames}
\title{
Write Reading Frames of Metagene to Files
}
\description{
Counting total number of reading frames for all metagene positions (positions of read start relative to cds start).
}
\usage{
getMetageneFrames(posProfile, bam_file_name)
}
\arguments{
  \item{posProfile}{
  A data frame with three columns for to_start,  to_end, and read_len of each reads in a Ba file.
}
  \item{bam_file_name}{
  Character vector, name of bam file from which the reads are scaned for metagene positions. Used for output file generation.
}
}
\details{
Two output files will be generated. One for metagene positions relative to cds start and one for metagene position relative to cds stop. Each file has three columns for total counts, frame, and metagene position.
}
\keyword{ methods }

