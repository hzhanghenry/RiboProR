\name{is.DNA.sequence}
\alias{is.DNA.sequence}
\title{
Check Out If A Sequence Is from DNA
}
\description{
Check out the given sequence is DNA by matching its based to A, T, G, and C.
}
\usage{
is.DNA.sequence(DNA_seq)
}
\arguments{
  \item{DNA_seq}{
  Character vector, s fragment of DNA sequence. No 'N' or white space allowed.
}
}
\value{
Logic, TRUE if the given sequence is from DNA Otherwise FALSE. 
}

\examples{
is.DNA.sequence("AGCTTAGGCCAAT") #	TRUE
is.DNA.sequence("AGCUUAGGCCAAU") #	FALSE
}
\keyword{ methods }
