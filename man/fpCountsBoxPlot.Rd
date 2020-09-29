\name{fpCountsBoxPlot}
\alias{fpCountsBoxPlot}
\title{
Counts Data Box Plot 
}
\description{
Make box plot with counts data fron DESeqDataSet Object
}
\usage{
fpCountsBoxPlot(dds, normalize=c("size", "fpm", "fpkm"), 
	label_area=2, main_text="Distribution of fp Counts")
}
\arguments{
  \item{dds}{
  A DESeqDataSet object
}
  \item{normalize}{
  Character vector, method of normalization either "size", "fpm", "fpkm", or NULL.
}
  \item{label_area}{
  Positive integer, height of area at botom of plot for sample labels. 
}
  \item{main_text}{
  Character vector, text for title of the plot.
}
}
\examples{
\dontrun{
data(dds_TE)
fpCountsBoxPlot(dds_TE)
}
}
\keyword{ methods }