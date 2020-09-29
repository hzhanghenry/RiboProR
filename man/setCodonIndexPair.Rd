\name{setCodonIndexPair}
\alias{setCodonIndexPair}
\title{
Pairing Start and Stop Codons in An ORF Range
}
\description{
Pair start and stop codons in an ORF ranges. For each start codon, only the  
first stop codon in its downstream can be used. If there is no stop codon  
found, the end of sequence is used (index 0).
}
\usage{
setCodonIndexPair(start_index, stop_index)
}
\arguments{
  \item{start_index}{
  Integer vector, start codon index, always greater than 0.
}
  \item{stop_index}{
  Integer vector, stop codon index,  could be 0 or integer(s) greater than 0.
}
}
\keyword{ methods }
