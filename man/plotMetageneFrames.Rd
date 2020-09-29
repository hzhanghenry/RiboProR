\name{plotMetageneFrames}
\alias{plotMetageneFrames}
\title{
Plot Frames of Metagene Regions
}
\description{
Make bar plot showing frame distribution at positions relative 
to CDS start and stop positions.  
}
\usage{
plotMetageneFrames(metagene_atStart=NULL, metagene_atStop=NULL,  
    min_5p=-20, max_5p=200, min_3p=-200, max_3p=20,  
    frame_colors=c("red", "green", "blue"), beside=TRUE)
}
\arguments{
  \item{metagene_atStart}{
Character vector, name of the file which contains frame name  
and frame counts for each metagene position relative to cds  
start position.  
}
  \item{metagene_atStop}{
Character vector, name of the file which contains frame name 
and frame counts for each metagene position relative to cds 
stop position.  
}
  \item{min_5p}{
  Integer, minimum distance to cds start positon.  
}
  \item{max_5p}{
  Integer, maximum distance to cds start positon.  
}
  \item{min_3p}{
  Integer, minimum distance to cds stop positon.  
}
  \item{max_3p}{
  Integer, maximum distance to cds stop positon.  
}
  \item{frame_colors}{
Character vector or R colors vector of length 3, plot colors 
for reading frames.  
}
  \item{beside}{
	Logic, if the columns for same group alongside each other. 
}
}
\examples{
\dontrun{
at_start <- system.file("data", "ribo_meta_start_frames.txt", 
					package="RiboProTools");
at_stopt <- system.file("data", "ribo_meta_stop_frames.txt", 
					package="RiboProTools");
plotMetageneFrames(at_start, at_stopt);
title("Metagene Frame Distribution");}
}
\keyword{ methods }
