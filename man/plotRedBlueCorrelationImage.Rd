\name{plotRedBlueCorrelationImage}
\alias{plotRedBlueCorrelationImage}
\title{
Plot Correlation Image with Red-Blue Colors
}
\description{
Plot pairwise correlation matrix with blue and red colors.
}
\usage{
plotRedBlueCorrelationImage(corr_data, image_name, 
	image_type="pdf", image_width=12)
}
\arguments{
  \item{corr_data}{
  Numeric matrix, pairwise correlation coefficients of samples.
}
  \item{image_name}{
  Character vector, output image file name.
}
  \item{image_type}{
  Character vector, output image format, either "pdf", "tiff", or "png".
}
  \item{image_width}{
  Positive integer, size of squared image.
}
}
\note{
There is no scree output and image will be saved to file.
}
\keyword{ methods }
