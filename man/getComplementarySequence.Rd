\name{getComplementarySequence}
\alias{getComplementarySequence}
\title{
Convert a fragment of sequence to its complementary contents
}
\description{
Get complementary sequence for a fragment of DNA sequence (no RNA)
}
\usage{
getComplementarySequence(seq_fragment)
}
\arguments{
  \item{seq_fragment}{
  Character vector, a fragment of DNA sequence.
}
}
\value{
  Character vector, the complementary sequence of the input with A -> T, G -> C, T -> A, and C -> G.
}
\keyword{ methods }
