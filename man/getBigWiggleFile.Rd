\name{getBigWiggleFile}
\alias{getBigWiggleFile}
\title{
Generate Coverage File in bigWig Format ()
}
\description{
Generate bigwig files from bam file. Only coverage of aligned reads are reported and no normalization performed. There will be two bigWig files saved to current directory, one for gorward strand and another for reverse strand.
}
\usage{
getBigWiggleFile(bam_file, transcripts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam_file}{
  Character vector, bam file name (and path)
}
  \item{transcripts}{
  GRange list for transcripts.
}
}
\author{
Henry Zhant
}
\examples{
	\dontrun{
	yeast_bed <- rtracklayer:import("yeast-all.bed");
	getBigWiggleFile("riboseq_t1.bam", yeast_bed)
	}
}
\keyword{ methods }
