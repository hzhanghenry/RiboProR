\name{getSequenceFromMutipleFastaFiles}
\alias{getSequenceFromMutipleFastaFiles}
\title{
Extract Genomic Sequences from Multiple Fasta Files
}
\description{
Read fasta files which are orgnaized by chromosome and extract relevant chromosome names and sequence.
}
\usage{
getSequenceFromMutipleFastaFiles(file_path, file_type)
}
\arguments{
  \item{file_path}{
character vector, path to the directory which holds the fasta files. 
}
  \item{file_type}{
Character vector, file extention. Valid types are "fa", "FA", "fasta", and "FASTA".
}
}
\value{
A data frame containing sequence data where rows are for each chromosome and columns are for chromosome names and sequences
}
\keyword{ methods }
