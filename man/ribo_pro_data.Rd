\name{ribo_pro_data}
\alias{ribo_pro_data}
\docType{data}
\title{
Sample Data for DESeq2 Analysis 
}
\description{
A numeric matrix of raw counts from RiboSeq and RNASeq Data including 5440 
genes and 6 samples for both RNASeq and RiboSeq. Row names are gene names.  
}
\usage{data("ribo_pro_data")}
\format{
  The format is: numeric matrix.
}
\source{
Unpublished data.
}

\examples{
\dontrun{
data(ribo_pro_data)
dim(ribo_pro_data)
head(ribo_pro_data)}
}
\keyword{datasets}
