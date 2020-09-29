\name{countAligns}
\alias{countAligns}
\title{
Count Aligned A-site for All Transcripts in Bed File
}
\description{
Scan a bam file and count aligned a-site for each transcript defined in bed file.
}
\usage{
countAligns(asiteOffsets, insets, bamfile, trxBed)
}
\arguments{
  \item{asiteOffsets}{
  A data frame of one column with read length as row names and a-sites for each read length in column.
}
  \item{insets}{
  List of integer with length of 4, insets in number of nucleotides for avoiding start and stop position.
}
  \item{bamfile}{
  Character vector, name of a bam file (and path).
}
  \item{trxBed}{
  GRange object for all transcripts.
}
}
\value{
List of list, number of aligned a-site in transcript, cds, 5'UTR, and 3'UTR region of each transcript.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }

