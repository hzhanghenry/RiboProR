\name{getSequence}
\alias{getSequence}
\title{
Extract Sequence From Fasta File
}
\description{
Get sequence from fasta file for a transcript. 
}
\usage{
getSequence(fafile, trx)
}
\arguments{
  \item{fafile}{
  Character vector, the fasta file name (and path).
}
  \item{trx}{
  GRanges object for a transcript.
}
}
\value{
Character vector, DNA sequence of the transcript.
}
\keyword{ methods }
