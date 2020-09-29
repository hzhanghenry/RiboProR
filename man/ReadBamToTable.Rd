\name{ReadBamToTable}
\alias{ReadBamToTable}
\title{
Scan A BAM File and Convert Results to A Data Frame
}
\description{
Read a bam file and convert the object (list of list) to a data frame
}
\usage{
ReadBamToTable(bam.file)
}
\arguments{
  \item{bam.file}{
  character vector, name of the bam file (and path).
}
}
\value{
A data frame with all conternts read from BAM file.
}
\references{
Morgan M, Pages H, Obenchain V and Hayden N (2016). Rsamtools: Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import. R  package version 1.26.1,  
 
http://bioconductor.org/packages/release/bioc/html/Rsamtools.html.
}
\keyword{ methods }
