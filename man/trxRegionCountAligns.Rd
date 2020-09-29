\name{trxRegionCountAligns}
\alias{trxRegionCountAligns}
\title{
Count Aligned A-site in A Transcripts by Regions
}
\description{
Check out if reads of a transcript has A site. If yes, get the total counts for transcript, cds, 5-UTR, and 3-UTR. Otherwise, set counts of the transcript to 0, and 0 or NA( if no cds defined) for cds, 5-UTR, and 3-UTR.
}
\usage{
trxRegionCountAligns(asiteOffsets, insets, bamfile, trx, cds)
}
\arguments{
  \item{asiteOffsets}{
  A data frame of one column with read length as row names and a-sites for each read length in column.
}
  \item{insets}{
  List of integer of length 4, insets in number of nucleotides to avoid start and stop cpdon.
}
  \item{bamfile}{
  Character vector, name of a bam file (and path).
}
  \item{trx}{
  GRange object for a transcript.
}
  \item{cds}{
  IRange object for cds in the transcript.
}
}
\value{
List of integer with length 4, number of aligned a-site in transcript, cds, 5-UTR, and 3-UTR regions.
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }
