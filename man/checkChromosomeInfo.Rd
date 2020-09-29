\name{checkChromosomeInfo}
\alias{checkChromosomeInfo}
\title{
Check Chromosome Information Between Fasta and BED Files
}
\description{
Check out if the chromsome names and length are same in the given two files.
}
\usage{
checkChromosomeInfo(DNA_seq, bed_info)
}
\arguments{
  \item{DNA_seq}{
  A data frame with rows for chromosome(s) and columns for chromosome name(s)   
  and sequence(s).
}
  \item{bed_info}{
A data frame with contents same as bigBed file.
}
}
\details{
The first file is a fasta file with chromosome(s) and sequences. 
The second file is an annotation file in BigBed format (12 columns).
This function will check if the number of chromosomes and chromosome names  
in the two files are same, and sequence length in fastq file is same as the  
chromosome length defined in bed file.
}
\value{
Logic, TRUE if the relevant information in two files are same. Otherwise return   
an error message.
}
\keyword{ methods }
