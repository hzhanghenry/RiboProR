\name{getCodonFromSequence}
\alias{getCodonFromSequence}
\title{
Split the Given mRNA Sequence to Codons
}
\description{
Covert the given mRNA sequence to codons starting the position of start_at.
}
\usage{
getCodonFromSequence(mRNA_seq, start_at)
}
\arguments{
  \item{mRNA_seq}{
  Character vector, mRNA seqwence (base A, U, G, C only).
}
  \item{start_at}{
  Positive integer, the start position to split the sequence. Must be 1, 2, or 3.
}
}
\value{
Character vector, codons starting from the start_at position in the sequence.
}
\examples{
getCodonFromSequence("GAUGAGCUAGGAC", start_at=2)
#	"AUG" "AGC" "UAG" "GAC"
}
\keyword{ methods }
