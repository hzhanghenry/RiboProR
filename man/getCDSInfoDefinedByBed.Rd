\name{getCDSInfoDefinedByBed}
\alias{getCDSInfoDefinedByBed}
\title{
Extract CDS Information and Sequence
}
\description{
Extract genomice sequence and other annotation items for CDS defined in bed file.
}
\usage{
getCDSInfoDefinedByBed(DNA_seq, bed_info)
}
\arguments{
  \item{DNA_seq}{
A data frame with rows for chromosome(s) and columns for chromosome name(s) and   
their sequence(s).
}
  \item{bed_info}{
A data frame with contents same as bigBed file (12 columns)
}
}
\value{
A data frame with 7 columns for:  
chromosome,  
start_pos,  
end_pos,  
strand,  
locus,  
sequence,  
type ("CDS")
}
\references{
https://genome.ucsc.edu/FAQ/FAQformat#format1.7
}
\keyword{ methods }
