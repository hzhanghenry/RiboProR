\name{covertToRNASequence}
\alias{covertToRNASequence}
\title{
Convert DNA Sequence to RNA
}
\description{
Convert DNA sequence to RNA sequence and reverse it if it is for reverse strand.
}
\usage{
covertToRNASequence(DNA_seq, strand)
}
\arguments{
  \item{DNA_seq}{
  Character vector, a fragment of DNA sequences from forward strand.
}
  \item{strand}{
  Character, either "+" or "-".
}
}
\value{
A character vector containing series of "A", "U", "G", "C".
}
\keyword{ methods }
