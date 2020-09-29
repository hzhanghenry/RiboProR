\name{filterOutAlignments}
\alias{filterOutAlignments}
\title{
Filter Aligments/Reads
}
\description{
Filter out alignments/reads by read length and genomic regions.
}
\usage{
filterOutAlignments(alignments, asite_table, annot_bed_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignments}{
  GRange list for alignments/reads in GRange list.
}
  \item{asite_table}{
  A data frame with one column for A-sites and rownames for qualified read length.
}
  \item{annot_bed_file}{
  Character vector, name (and path) of annotation file in bed format.
}
}
\value{
 GRange list with qualified alignments/reads.
}
\keyword{ methods } 
