\name{quickDESeq2Test}
\alias{quickDESeq2Test}
\title{
Perform A Quick DESeq Analysis for Ribosomal Footprints Profiling 
}
\description{
A simple way to perform DESeq2 analysis with ribosomal footprints profiling data and gene annotation information.
}
\usage{
quickDESeq2Test(count_table, Ribo_wildtype, Ribo_mutants, 
		mRNA_wildtype, mRNA_mutants, control_name, mutant_name, 
		mRNA_level, Ribo_level, annotation_file)
}
\arguments{
  \item{count_table}{
  Numeric matrix with RNASeq reads counts and Riboseq footprints counts for both control and treamtment/mutant samples.
}
  \item{Ribo_wildtype}{
  Positive integer vectors, column numbers in count matrix for raw counts of RiboSeq wildtype samples.
}
  \item{Ribo_mutants}{
  Positive integer vectors, column numbers in count matrix for raw counts of RiboSeq treatment/mutant samples.
}
  \item{mRNA_wildtype}{
  Positive integer vectors, column numbers in count matrix for raw counts of RNASeq wildtype samples.
}
  \item{mRNA_mutants}{
  Positive integer vectors, column numbers in count matrix for raw counts of RNASeq treatment/mutant samples.
}
  \item{control_name}{
  Character vector, name of control group.
}
  \item{mutant_name}{
  Character vector, names of treatment/mutant group.
}
  \item{mRNA_level}{
  Numeric, threshold to filter out the matrix based on row means of RNASeq samples.
}
  \item{Ribo_level}{
  Numeric, threshold to filter out the matrix based on row means of RiboSeq samples.
}
  \item{annotation_file}{
  Character vector, name of gene annotation file
}
}
\value{
None. All results will be saved to file.
}
\examples{
\dontrun{
data("ribo_pro_data.RData")
annotation_file <- system.file("data", "yeast_gene_descriptions.txt", 
		package="RiboProTools");
quickDESeq2Test(ribo_pro_data, Ribo_wildtype=c(1:3), Ribo_mutants=c(7:9), 
		mRNA_wildtype=c(4:6), mRNA_mutants=c(10:12), 
		control_name="ribo_WT", mutant_name="ribo_MT", 
		mRNA_level=10, Ribo_level=1, annotation_file)
}
}
\keyword{ methods }

