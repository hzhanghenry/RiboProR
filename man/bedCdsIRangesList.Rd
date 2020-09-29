\name{bedCdsIRangesList}
\alias{bedCdsIRangesList}
\title{
Get IRanges List for All CDS
}
\description{
For each transcript in a BED file and the associated transcript GRanges as computed by bedGRangesList, find the transcript-local CDS as per transcriptCdsIRanges and collect these as an IRangesList.
}
\usage{
bedCdsIRangesList(bed, bedgrl)
}
\arguments{
  \item{bed}{
  GRanges list imported from bed file.
}
  \item{bedgrl}{
  GRanges list of exons for all transcripts.
}
}
\value{
IRanges list for each exon (CDS).
}
\references{
Original code from Nicholas T. Ingolia, et al. (2014). Ribosome Profiling Reveals Pervasive Translation Outside of Annotated Protein-Coding Genes. Cell Reports 8, 1365-1379.
}
\keyword{ methods }

