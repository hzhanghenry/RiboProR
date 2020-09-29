\name{writeChromosomeSizesToFile}
\alias{writeChromosomeSizesToFile}
\alias{readChromosomeSizesFromFile}
\title{
Write Chromosome Sizes to File
}
\description{
Read in chromosome names and sizes from or Write chromosome names and sizes to a tab-delimited text file.
}
\usage{
writeChromosomeSizesToFile(chromSizes, file_name)
readChromosomeSizesFromFile(file_name)
}
\arguments{
  \item{chromSizes}{
  A data frame with 2 columns for chromosome names and lengths.
}
  \item{file_name}{
  Character vector, name of file (and path) for read from or write to.
}
}
\value{
readChromosomeSizesFromFile() returns a data frame with 2 columns for chromosome names and lengths.
}
\keyword{ methods }

