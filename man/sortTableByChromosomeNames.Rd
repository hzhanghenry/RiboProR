\name{sortTableByChromosomeNames}
\alias{sortTableByChromosomeNames}
\title{
Sort A Table by Chromosome Names
}
\description{
Sort a table by chromosome names (either Arabic or Roman numbers)
}
\usage{
sortTableByChromosomeNames(chrom_info, name_col=1, type="digit")
}
\arguments{
  \item{chrom_info}{
  A data frame or matrix with one column for chromosome names.
}
  \item{name_col}{
  Positive integer, number of the column for chromosome names.
}
  \item{type}{
  Character vector, type of chromosome numbers, either "digit" or "roman".
}
}
\value{
Data frame or matrix same as the input but sorted by chromosome names.
}
\keyword{ methods }
