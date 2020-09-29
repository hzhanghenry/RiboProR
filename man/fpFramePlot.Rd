\name{fpFramePlot}
\alias{fpFramePlot}
\title{
Bar Plot of Reading Frames Data
}
\description{
Plot reading frames distributions of ribosomal foot prints (bar plot).
}
\usage{
fpFramePlot(frame_file, by_column=TRUE, fp_column=c(2:4), read_len=c(25:30), 
	legend_pos="topleft", title_text="Ribosomal fp Frame Distribution")
}
\arguments{
  \item{frame_file}{
  Character vector, name of the file (and path) with reading frame counts for each read length.
}
  \item{by_column}{
  Logic, if the frame data arraned by column (row is for each read length). If TRUE, fp_column must be defined.
}
  \item{fp_column}{
  Positive integer, column number for reading frame counts.
}
  \item{read_len}{
  Positive integer vector, length of reads (alignments).
}
  \item{legend_pos}{
  Character vector, location of legend, one of "topleft", "center", "topright".
 }
  \item{title_text}{
  Character vector, text for title.
}
}
\examples{
\dontrun{
frame_file <- system.file("extdata", "RiboProTools_ribo.frame_len.txt",  
	package = "RiboProTools")
fpFramePlot(frame_file)
}
}
\keyword{ methods }

