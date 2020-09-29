\name{getCountMatrixFromTable}
\alias{getCountMatrixFromTable}
\title{
Extract Raw Counts from A Matrix or Data Frame
}
\description{
Generate a matrix from a matrix or data frame with defined columns.
}
\usage{
getCountMatrixFromTable(count_table=NULL, ribo_control=NULL, 
	ribo_treatment=NULL, mRNA_control=NULL, mRNA_treatment=NULL)
}

\arguments{
  \item{count_table}{
  A data frame or matrix with raw counts from both Riboseq and RNAseq data for same samples. Row names must be gene names or gene IDs.
}
  \item{ribo_control}{
  Positive integer vector, columns in input table for control samples of Riboseq data.
}
  \item{ribo_treatment}{
  Positive integer vector, columns in input table for treatment or mutant samples of Riboseq data.
}
  \item{mRNA_control}{
  Positive integer vector, columns in input table for control samples of RNASeq data.
}
  \item{mRNA_treatment}{
  Positive integer vector, columns in input table for treatment or mutant samples of RNASeq data.
}
}
\value{
A matrix with columns in the order of mRNA control samples, mRNA treatment samples, ribosomal control samples, and ribosomal treatment samples
}
\examples{
\dontrun{
data("ribo_pro_data.RData")
getCountMatrixFromTable(count_table=ribo_pro_data, ribo_control=7:9, 
	ribo_treatment=10:12, mRNA_control=1:3, mRNA_treatment=4:6)
}
}
\keyword{ methods }
