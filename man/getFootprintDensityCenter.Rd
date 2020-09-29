\name{getFootprintDensityCenter}
\alias{getFootprintDensityCenter}
\title{
Main Function to Calculate Center of Ribosomal Footprint Densities 
}
\description{
Calculate centers of ribosomal footprint density for each transcript in a bam file .
}
\usage{
getFootprintDensityCenter(a_sites, transcripts, bam_file, weight = FALSE)
}
\arguments{
  \item{a_sites}{
Data frame with 1 column for a-site. Row names are read length.
}
  \item{transcripts}{
GRanges list of all transcripts.
}
  \item{bam_file}{
Character vector, name (and path) of a bam file.
}
  \item{weight}{
Logic, if the position with zero count should be filtered out.
}
}
\value{
A data frame with one column of positive integers for density center of all transcripts. Row names of the data frame are gene ID/names.
}
\keyword{ methods }
