\name{getCountsData}
\alias{getCountsData}
\title{
Get Counts Data from DESeqDataSet Object
}
\description{
Extract raw counts data or normalized counts from DESeqDataSet object.
}
\usage{
getCountsData(dds, normalize=c("size", "fpm", "fpkm"))
}
\arguments{
  \item{dds}{
  A DESeqDataSet object.
}
  \item{normalize}{
  Character vector, method of normalization either "size", "fpm", "fpkm", or NULL.
}
}
\value{
A matrix holds raw or normalize counts data.
}
\examples{
\dontrun{
data("dds_TE")
fpm <- getCountsData(dds_TE, normalize="fpm")
}
}
\keyword{ methods }

