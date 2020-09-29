\name{plotMetageneCounts}
\alias{plotMetageneCounts}
\title{
Plot Metagene Counts Distribution
}
\description{
Make line plot showing metage counts distribution around positions relative 
to CDS start and stop positions.  
}
\usage{
plotMetageneCounts(count_file, from_position=-50, 
	to_position=50, codon="Start", x_interval=10)
}
\arguments{
  \item{count_file}{
Character vector, name of the file which contains frame name  
and frame counts for each metagene position relative to cds  
start position.  
}
  \item{from_position}{
  Integer, leftmost position relative to cds start or stop position.  
}
  \item{to_position}{
  Integer, rightmost position relative to cds start or stop position.  
}
  \item{codon}{
Character vector, either 'Start' or 'Stop'.  
}
  \item{x_interval}{
	Integer, length between tick-marks for x-axis. 
}
}
\examples{
\dontrun{
at_start <- system.file("data", "ribo_meta_start_frames.txt", 
plotMetageneCounts(at_start);
title("Metagene Counts Distribution");}
}
\keyword{ methods }
