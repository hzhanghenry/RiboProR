\name{getPrettyLabels}
\alias{getPrettyLabels}
\title{
Format Axis Labels
}
\description{
Modifiy the axis labels for metagene frame plot.  
}
\usage{
getPrettyLabels(min_val, max_val)
}
\arguments{
  \item{min_val}{
  Integer, minimum value of the axis label.  
}
  \item{max_val}{
  Integer, maximum value of the axis label.  
}
}
\value{
Integer vector.  
}
\references{
http://bioconductor.org/packages/riboSeqR/
}
\examples{
label_5p <- getPrettyLabels(-50, 200);
}
\keyword{ methods }
