\name{regionCountFrame}
\alias{regionCountFrame}
\title{
Convert Aligned A-site of All Transcripts to A Data Frame
}
\description{
Convert the aligned a-site of all transcripts from a list of list to a data frame.
}
\usage{
regionCountFrame(counts)
}
\arguments{
  \item{counts}{
List of list, number of aligned a-site in transcript, cds, 5'UTR, and 3'UTR region of each transcript.
}
}
\value{
Data frame with 4 columns, total number of aligned a-site in transcript, cds, 5'UTR, and 3'UTR region of each transcript.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
