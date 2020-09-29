\name{plotReplicates}
\alias{plotReplicates}
\title{
Plot Raw or Normalized Counts of Sample Replicates
}
\description{
Make scatter plot with count data of two replicate samples in DESeqDataSet.
}
\usage{
plotReplicates(dds, normalize=c("size", "fpm", "fpkm"), 
	replicates=c(1, 2), show.cor=TRUE, is.log2=FALSE, 
	point_color="grey", line_color="red", 
	x_label="replicate 1", y_label="replicate 2", 
	main_text="Distribution of fp Counts")
}
\arguments{
  \item{dds}{
  A DESeqDataSet object.
}
  \item{normalize}{
  Character vector, method of normalization either "size", "fpm", "fpkm", or NULL.
}
  \item{replicates}{
  Positive integer vector of length 2. The columns in counts data to be plotted.
}
  \item{show.cor}{
  Logic, if plot corelation coefficiency and the regression line.
}
  \item{is.log2}{
  Logic, is the data log2 transformed.
}
  \item{point_color}{
  Character vector of R colors names, color for the point background.
}
  \item{line_color}{
  Character vector, colors for regression line.
}
  \item{x_label}{
  Character vector, text for x axis label
}
  \item{y_label}{
  Character vector, text for y axis label
}
  \item{main_text}{
  Character vector, text for title of the plot
}
}
\examples{
\dontrun{
data("dds_TE")
plotReplicates(dds_TE, replicates=1:2)
}
}
\keyword{ methods }
