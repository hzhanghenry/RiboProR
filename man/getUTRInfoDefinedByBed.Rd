\name{getUTRInfoDefinedByBed}
\alias{getUTRInfoDefinedByBed}
\title{
Extract UTR Information From Both of Fasta and BED Files
}
\description{
Extract annotation information from bed file and sequence from fasta file.
}
\usage{
getUTRInfoDefinedByBed(DNA_seq, bed_info)
}
\arguments{
  \item{DNA_seq}{
A data frame with rows for chromosome(s) and columns for chromosome name(s) and sequence.
}
  \item{bed_info}{
A data frame with contents same as bigBed file.
}
}
\value{
A data frame with 7 columns for:

chromosome:	chromosome name of each UTR 
start_pos:	start position of each UTR 
end_pos:	end position of each UTR 
strand: 	strand of each UTR
locus: 		gene name of each UTR
sequence:	DNA sequence of each UTR
type:		type of each fragment (5UTR or 3UTR)

All positions are forward strand based.
}
\references{
https://genome.ucsc.edu/FAQ/FAQformat#format1.7
}
\keyword{ methods }
