\name{getCodonIndex}
\alias{getCodonIndex}
\title{
Get the Index of A Codon in Codon List
}
\description{
Find which codon(s) in the codon list matches to the target codon.
}
\usage{
getCodonIndex(codon_list, target_codon)
}
\arguments{
  \item{codon_list}{
Character vector, a serial codons from a RNA sequences.
}
  \item{target_codon}{
Character vector, one or more start or stop codon(s).
}
}
\value{
Positive integer vector, index of the codon(s) which matched to target codon.
}
\keyword{ methods }
