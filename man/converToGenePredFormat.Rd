\name{converToGenePredFormat}
\alias{converToGenePredFormat}
\title{
Covert ORF Information to GenePred Format
}
\description{
Covert ORF table in tab-delimited format to GenePrep format.
}
\usage{
converToGenePredFormat(ORF_file, is.one.based = TRUE, min_codon = 0)
}
\arguments{
  \item{ORF_file}{
Character vector, name (and path) of tab-delimited file for ORF information.
}
  \item{is.one.based}{
Logic, if position is 1-based in input file. Default TRUE.
}
  \item{min_codon}{
Integer. 0 or bigger. Number of codon(s) between start and stop codons.
}
}
\details{
The output has following columns (one column more than the ucsc gebe prediction  
table formt:
gene_name  
orf_ID  
chromosomes  
strand  
tx_start  
tx_end  
cds_start  
cds_end  
num_exon  
exon_starts  
exon_end  
}
\value{
Data frame with 11 columns in GenePred format.
}
\references{
https://genome.ucsc.edu/FAQ/FAQformat.html#format9
}
\keyword{ methods } 