\name{getPlotColors}
\alias{getPlotColors}
\alias{getDefaultColors}
\title{
Get R Colors for Plot
}
\description{
Generate a default color list of length 24 or rainbow colors with length more than 24.
}
\usage{
getPlotColors(num_colors)
getDefaultColors()
}
\arguments{
  \item{num_colors}{
  Positive integer, total number of colors.
}
}
\value{
getPlotColors(num_colors) returns a vector of R colors with length of num_colors.  

getDefaultColors() returns a vector of predifined R colors with length of 24.
}
\keyword{ methods }
