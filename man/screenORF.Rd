\name{screenORF}
\alias{screenORF}
\title{
Extract ORF Information from A Sequence Fragment
}
\description{
Screen a sequence fragment to get start codon, context from -3  to +4 base  
positions of start codon, length of ORF, distance from the cap, and distance  
to the uAUG in a list.
}
\usage{
screenORF(seq_info, codon_list, codon_pair, start_at, mRNA_seq)
}
\arguments{
  \item{seq_info}{
One row from a data frame, sequence info including of chromosome, start_pos,  
end_pos, strand, locus, sequence, and type.
}
  \item{codon_list}{
Character vector, codons converted from UTR sequence.
}
  \item{codon_pair}{
Matrix of column 2 for paired start and stop index.
}
  \item{start_at}{
Integer in 1, 2, 3, base where start to read (frame).
}
  \item{mRNA_seq}{
Character vector, a fragment of mRNA sequence.
}
}
\details{
Screen a mRNA sequence to get information including of:  
gene_name  
chromosome  
seq_start  
seq_stop  
strand  
orf_id  
orf_start  
orf_stop 
start_codon 
context   
orf_length  
dis_from_cap  
dis_to_MUG  
}
\value{
A data frame with columns for information above.
}
\keyword{ methods }
