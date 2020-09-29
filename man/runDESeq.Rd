\name{runDESeq}
\alias{runDESeq}
\title{
Run DESeq on A DESeqDataSet
}
\description{
Run DESeq() or other related functions based on the contents of the DESeqDataSet.
}
\usage{
runDESeq(deseq_dataset, has.SizeFactors=FALSE, 
	reset.design=FALSE, fit_type="parametric")
}
\arguments{
  \item{deseq_dataset}{
  A DESeqDataSet with design model of ~ genotype + condition + genotype:condition, on which DESeq() or other related functions will be called.
}
  \item{has.SizeFactors}{
  Logic, if the deseq_dataset has size factor already or not.
}
  \item{reset.design}{
  Logic, if need reset design model.
}
  \item{fit_type}{
  Character vector, fit type used by DESeq() or estimateDispersions().
}
}
\details{
This function will apply DESeq() or other related functions based on the contents of the DESeqDataSet object. If the DESeqDataSet has no sizeFactors calculated, it will call DESeq() otherwise it will call estimateDispersions() and nbinomWaldTest(). By default, the design model is ~ genotype + condition + genotype:condition, which is for translational efficiency test. For tranascription test, it must be changed to ~group. This will be done if set reset.design to TRUE.
}
\value{
A DESeqDataSet, on which DESeq() or relative functions has been called.
}
\references{
https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf
}
\keyword{ methods }
