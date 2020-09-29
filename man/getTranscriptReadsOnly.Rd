\name{getTranscriptReadsOnly}
\alias{getTranscriptReadsOnly}
\title{
Filter Alignments/Reads in GRange List by Genomic Regions
}
\description{
Filter out alignments/reads to keep that overlapped with genomic regions defined in annotation file.
}
\usage{
getTranscriptReadsOnly(alignments, annot_bed_file)
}
\arguments{
  \item{alignments}{
  GRange object of all alignments to be filtered out.
}
  \item{annot_bed_file}{
  Character vector, name (and path) of annotation file in bed format.
}
}
\value{
GRange object of filtered alignments.
}
\keyword{ methods }% __ONLY ONE__ keyword per line
