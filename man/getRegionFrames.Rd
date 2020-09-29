\name{getRegionFrames}
\alias{getRegionFrames}
\title{
Generate Metagene Frame Table for A Specific Region
}
\description{
Extract a subset from full metagene frame table and convert them to a three row data frame.  
}
\usage{
getRegionFrames(frame_file, from, to)
}
\arguments{
  \item{frame_file}{
  Character vecor, name (and path) of the file continning frames 
  at each metagene position.  
}
  \item{from}{
  Integer, the start of a sub region of the metagene.  
}
  \item{to}{
  Integer, the end of a sub region of the metagene.  
}
}
\value{
A data frame with 3 rows for counts of frame0, frame1, and frame2 at each metagene position.  
}
\examples{
\dontrun{
frame_file <- system.file("extdata", "frame_at_start.txt", 
	package="RiboProTools")
frame_table <- getRegionFrames(frame_file, -50, 200);
}
}
\keyword{ methods }
