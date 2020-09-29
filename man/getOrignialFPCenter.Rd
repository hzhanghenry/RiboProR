\name{getOriginalFPCenter}
\alias{getOriginalFPCenter}
\title{
Calculate Center of Ribosomal Footprint Density Without Filtering
}
\description{
Calculate center of ribosomal footprint density for a transcript without filtering.
}
\usage{
getOriginalFPCenter(asiteProfile, trx_length)
}
\arguments{
  \item{asiteProfile}{
Vector of positive integers, count of asite aligned to the transcript at each nt position.
}
  \item{trx_length}{
Positive integer, length of the transcript.
}
}
\value{
Positive float number, ratio of length of the 5' end of transcript that has half of total counts divided by transcript length
}
\keyword{ methods }
