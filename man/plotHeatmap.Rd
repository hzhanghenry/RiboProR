\name{plotHeatmap}
\alias{plotHeatmap}
\title{
Heatmap Plot with Data Matrix
}
\description{
Make heatmap plot with data matrix and heatmap.2() in gplots package.
}
\usage{
plotHeatmap(plot_value, sample_name, gene_name, 
	image_type="pdf", image_width = 8, 
	is.log2=FALSE, scale_by="row")
}
\arguments{
  \item{plot_value}{
  Numeric matrix, values for heatmap plot.
}
  \item{sample_name}{
  Character vector, column (sample) labels.
}
  \item{gene_name}{
  Character vector, row (gene) labels.
}
  \item{image_type}{
  Character vector, output image format, one of "pdf", "tiff", and "png".
}
  \item{image_width}{
  Positive integer, width of output image in inches.
}
  \item{is.log2}{
  Logic, is the data log2 transformed.
}
  \item{scale_by}{
  Character vector, how the data is scaled, either "row" or "column"
}
}
\keyword{ methods }
