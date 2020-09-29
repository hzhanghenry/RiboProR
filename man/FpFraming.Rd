\name{FpFraming}
\alias{FpFraming}
\title{
Get Framing Information for All Reads in A Bam File
}
\description{
Perform framing analysis for all reads in a bam file and generate profiles for reading frames, metagenes, and a-sites. 
}
\usage{
FpFraming(bam_file, bed_file, parameters=NULL, validate.parameters=TRUE, save.file=TRUE)
}
\arguments{
  \item{bam_file}{
  Character vector, name of bam file (and path).
}
  \item{bed_file}{
  Character vector, name of bed file (and path) for gene annotation.
}
  \item{parameters}{
  List of integers for FpFraming parameters.
}
  \item{validate.parameters}{
  Logic, if TRUE, non-default parameters will be validated.
}
  \item{save.file}{
  Logic, if true, all tables (total of 6) will be saved to files.
}
}
\details{
This is the main function for framing analysis. It scans a bam file and extract information for metagene, read frames, and find a-sites for reads of each qualified length. Outputs will be saved as tab-delimited text files.
}
\value{
A data frame containing a-sites for fragments with selected length.
}
\references{
https://github.com/ingolia-lab/RiboSeq
}
\keyword{ methods }
