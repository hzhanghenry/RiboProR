\name{plotCorrelationHeatmap}
\alias{plotCorrelationHeatmap}
\title{
Plot Correlation Heastmap for samples
}
\description{
Make a correlation image for samples in DESeqDataSet with blue and red colors.
}
\usage{
plotCorrelationHeatmap(dds, normalize=c("size", "fpm", "fpkm"), 
	image_name, image_type="pdf", image_size=12)
}
\arguments{
  \item{dds}{
  A DESeqDataSet object.
}
  \item{normalize}{
  Character vector, method of normalization either "size", "fpm", "fpkm", or NULL.
}
  \item{image_name}{
  Character vector, output image file name.
}
  \item{image_type}{
  Character vector, output image format, either "pdf", "tiff", or "png".
}
  \item{image_size}{
  Positive integer, image height.
}
}
\examples{
\dontrun{
data("dds_TE")
plotCorrelationHeatmap(dds_TE, image_name="corr.image.pdf")
}
}
\keyword{ methods }
