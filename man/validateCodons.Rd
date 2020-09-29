\name{validateCodons}
\alias{validateCodons}
\title{
Check Out If the Given Codons Are Qualified Ones
}
\description{
Check out if each codon has correct length and correct bases.
}
\usage{
validateCodons(codons)
}
\arguments{
  \item{codons}{
character vector which holds one or more codons.
}
}
\details{
Each codon must have exactly 3 bases from A, U, G, and C.
}
\value{
None. Error message will be generated if any codon is a invalid one.
}
\keyword{methods }
