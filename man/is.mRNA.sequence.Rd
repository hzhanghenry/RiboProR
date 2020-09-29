\name{is.mRNA.sequence}
\alias{is.mRNA.sequence}
\title{
Check Out If A Sequence Is from mRNA
}
\description{
Check out if the given sequence is from mRNA, i.e., its bases must be A, U, G, or C.
}
\usage{
is.mRNA.sequence(mRNA_seq)
}
\arguments{
  \item{mRNA_seq}{
  Character vector, s fragment of mRNA sequence. No 'N' or white space allowed.
}
}
\value{
Logic, TRUE if the given sequence is from mRNA Otherwise FALSE. 
}
\examples{
is.mRNA.sequence("AUGCGAAUGGCC")	#	TRUE
is.mRNA.sequence("ATGCGAATGGCC")	#	FALSE
}
\keyword{ methods }
