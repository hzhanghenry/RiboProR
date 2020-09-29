\name{getSequenceFromOneFastaFile}
\alias{getSequenceFromOneFastaFile}
\title{
Extract Chromosome Names and Sequences from Fasta File
}
\description{
Read fasta file which include multiple chromosomes and relevant sequences. Simply put them in data frame for easy use.
}
\usage{
getSequenceFromOneFastaFile(file_name)
}
\arguments{
  \item{file_name}{
Character vector, name (and path) of a fasta file.
}
}
\value{
A data frame wiht row(s) for chromosome(s) and columns for chromosome names and sequences.
}
\keyword{ methods }
