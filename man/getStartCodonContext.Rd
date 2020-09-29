\name{getStartCodonContext}
\alias{getStartCodonContext}
\title{
Extract Contex for A Start Codon
}
\description{
Build the context (a short sequence fragment from the third nucleotide before  
start codon and the one nucleotide next to start codon). There will be total   
of 7 nucleotides, e.g., GACAUGG, AUCAUGC).
}
\usage{
getStartCodonContext(codon_list, start_index, mRNA_seq, start_at)
}
\arguments{
  \item{codon_list}{
Character vector, codons from the sequence.
}
  \item{start_index}{
Integer, index of start codon. Always greater than 0.
}
  \item{mRNA_seq}{
Character vector, mRNA sequence converted from UTR/CDS sequence
}
  \item{start_at}{
Integer, either 1, 2 or 3.
}
}
\value{
Character vector, nucleotides around start codon.  
}
\keyword{ methods }
