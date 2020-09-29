\name{extractEfficiencyChange}
\alias{extractEfficiencyChange}
\title{
Extract Results from DESeqDataSet Object For Transclational Efficiency Change.
}
\description{
Extract translational efficiency (TE) and transclational efficiency change (TEC) from DESeqDataSet object after calling of DESeq().  The output will contain all columns from results() and mcols().
}
\usage{
extractEfficiencyChange(dds_object, control_name, mutant_name, 
	efficiency_type="TE", meta_cols=c(1:3))
}
\arguments{
  \item{dds_object}{
  A DESeq object, on which DESeq() has already been called. The design model must be ~ genotype + condition + genotype:condition.
}
  \item{control_name}{
  Character vectors, name of control group, used for generation of output file name.
}
  \item{mutant_name}{
  Character vectors, name of mutant group, used for generation of output file name.
}
  \item{efficiency_type}{
  Character vectors, name of efficiency. Use "TE" for translational efficiency analysis or "RRO" for relative ribosomal occupancy analysis.
}
  \item{meta_cols}{
  Positive integers, columns of meta-data to be attached to output. Set to NULL if all columns are attached or 0 for no meta-data.
}
}
\value{
DESeqDataSet object, same for input. All results extracted are saved to files.
}
\references{
https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
}
\keyword{ methods }
