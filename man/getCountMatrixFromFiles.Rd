\name{getCountMatrixFromFiles}
\alias{getCountMatrixFromFiles}
\title{
Read Raw Counts from Text Files and Generate Counts Matrix
}
\description{
Read multiple tab-delimited text files and extract raw counts column, then merge all raw counts as one matrix with rows for genes and columns for samples.
}
\usage{
getCountMatrixFromFiles(directory_name, file_name_pattern, count_column, 
	rowname_column = 0, has.header = TRUE)
}
\arguments{
  \item{directory_name}{
  character vector, name of directory where raw counts files are stored.
}
  \item{file_name_pattern}{
  character vector, common pattern in all raw counts files. This pattern should not be used by other files in same directory.
}
  \item{count_column}{
  Positive integer, number of column for raw counts in the file.
}
  \item{rowname_column}{
   Positive integer, the column in raw count file for row names. Set to 0 for raw count files saved by write.table() with row.names=TRUE.
}
  \item{has.header}{
  Logical, if the raw count files have column headers.
}
}
\value{
A numeric matrix with rows for genes and columns for samples. File names are used as column headers in this matrix.
}
\note{
All raw counts files must be generated with same software or have same number of rows and columns as well as same row and column orders).
}
\examples{
directory_name <- "raw_count_files"
file_name_pattern <- "raw_counts.txt"
count_column <- 2
rowname_column <- 0
has.header <- TRUE

\dontrun{getCountMatrixFromFiles(directory_name, file_name_pattern, 
	count_column, rowname_column, has.header)}
}
\keyword{ methods }
