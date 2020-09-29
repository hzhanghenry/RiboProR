\name{extractTranscriptionChange}
\alias{extractTranscriptionChange}
\title{
Extract Results from DESeqDataSet Object For Transcription Change
}
\description{
Extract transcription change between Riboseq samples and between RNASeq samples from DESeqDataSet object after running of DESeq().   The output will contain all columns from results() and mcols().
}
\usage{
extractTranscriptionChange(dds_object, control_name, mutant_name, 
		change_one="mRNA", change_two="Ribo", meta_cols=c(1:3))
}
\arguments{
  \item{dds_object}{
  A DESeq object, on which DESeq() has already been called. The design model must use a grouping variable.
}
  \item{control_name}{
  Character vectors, name of control group, used for generation of output file name.
}
  \item{mutant_name}{
  Character vectors, name of mutant group, used for generation of output file name.
}
  \item{change_one}{
  Character vector, name of changes between groups of RNASeq samples.
}
  \item{change_two}{
  Character vector, name of changes between groups of RiboSeq samples.
}
  \item{meta_cols}{
  Positive integers, columns of meta-data to be attached to output. Positive integers, columns of meta-data to be attached to output. Set to NULL if all columns are attached or 0 for no meta-data.
}
}
\value{
A DESeq object, same as the input.
}
\references{
https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
}
\keyword{ methods }

