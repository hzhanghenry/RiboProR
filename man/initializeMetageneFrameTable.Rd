\name{initializeMetageneFrameTable}
\alias{initializeMetageneFrameTable}
\title{
Initialize A Metagene Frame Table
}
\description{
Initialized a matirx with 3 columns for counts, frames, and unique metagene positions.
}
\usage{
initializeMetageneFrameTable(position_set)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{position_set}{
  Integer vector, unique position set of all metagene positions (position relative to cds start or cds stop positions)
}
}
\value{
A data frame of 3 columns for number of frames, frame names, and metagene positions.
}
\references{
https://github.com/ingolia-lab/RiboSeq
}
\examples{
\dontrun{
positions <- -100:100;
frames_table <- initializeMetageneFrameTable(positions);
head(frames_table);
}
}
\keyword{ methods }
