\name{getFilteredFPCenter}
\alias{getFilteredFPCenter}
\title{
Calculate Center of Filtered Ribosomal Footprints for A Transcript
}
\description{
Calculate center of ribosomal footprints for a transcript by removing positions where have no a-site aligned.
}
\usage{
getFilteredFPCenter(asiteProfile)
}
\arguments{
  \item{asiteProfile}{
Vector of positive integer, counts of a-site on each nucleotide position of a transcript
}
}
\value{
Positive float number, length of the 5' end of transcript that has half of total counts divided by transcript length
}
\keyword{ methods }

