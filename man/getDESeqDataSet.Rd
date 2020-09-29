\name{getDESeqDataSet}
\alias{getDESeqDataSet}
\title{
Initialize A DESeqDataSet object
}
\description{
Generate a DESeqDataSet object with a counts matrix that contains control and treatment/mutant samples from both RiboSeq and RNASeq data.
}
\usage{
getDESeqDataSet(count_matrix, num_Ribo_wildtype, num_Ribo_mutant, 
	num_mRNA_wildtype, num_mRNA_mutant, annotation_info=NULL)
}
\arguments{
  \item{count_matrix}{
  A matrix with columns in the order of: mRNA control samples, and mRNA treatment samples, ribosomal control samples, and ribosomal treatment samples.
}
  \item{num_Ribo_wildtype}{
  Positive integer, total number of wildtype/control samples with RiboSeq fp counts.
}
  \item{num_Ribo_mutant}{
  Positive integer, total number of mutant/treament samples with RiboSeq fp counts.
}
  \item{num_mRNA_wildtype}{
  Positive integer, total number of wildtype/control samples with RNASeq fp counts.
}
  \item{num_mRNA_mutant}{
  Positive integer, total number of mutant/treament samples with RNASeq fp counts.
}
  \item{annotation_info}{
  Data frame with columns for gene ID, gene name, and description.
}
}
\details{
This function set up a DESeqDataSet using DESeqDataSetFromMatrix(). The countData will be the count_matrix, colData will have both condition ("mRNA and "Ribo") and genotype("wildtype" and "mutant"), and design will be ~ genotype + condition + genotype:condition. The annotation_info will be added into mcols, if provided.  
}
\value{
A DESeqDataSet object. 
}
\references{
https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
}
\keyword{ methods }

