\name{getORFPositions}
\alias{getORFPositions}
\title{
Calculate Relevant Genomic Positions for a ORF
}
\description{
Calculate genomic positions, for a ORF, includeing of start, stop, distance  
to cap (UTR start), and distance to main AUG codon.
}
\usage{
getORFPositions(seq_info, start_index, stop_index, start_at, include_stop_codon)
}
\arguments{
  \item{seq_info}{
Information of UTR (one row from data frame) 
}
  \item{start_index}{
Integer, start codon index in an ORF.
}
  \item{stop_index}{
Integer, stop codon index in an ORF.
}
  \item{start_at}{
Integer, start point in sequence to read codons (frame)
}
  \item{include_stop_codon}{
Logic, if include stop codon in outputs.
}
}
\value{
Positive integer vector for start and stop positions, distance from cap and   
to main AUG codon.
}
\keyword{ methods } 
