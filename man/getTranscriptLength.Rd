\name{getTranscriptLength}
\alias{getTranscriptLength}
\title{
Calculate Length of Transcripts
}
\description{
Calculate transcript length for exon only (remove introns, if any).
}
\usage{
getTranscriptLength(transcripts)
}
\arguments{
  \item{transcripts}{
GRange list with transcripts
}
}
\value{
A data frame of one column for transcript length and with gene names as row names.
}
\keyword{ methods }
