\name{getUTRInfoFromFastaFile}
\alias{getUTRInfoFromFastaFile}
\title{
Extract UTR Information from Fasta File
}
\description{
Process fasta file for UTR to get annotation items and sequence for each UTR,  
e.g., 5UTR.fa downloaded from https://www.pombase.org/downloads/utr
}
\usage{
getUTRInfoFromFastaFile(fasta_file, header_name)
}
\arguments{
  \item{fasta_file}{
  Character vector, name (and path) of fasta file.
}
  \item{header_name}{
  Character vector, definition of each header field.
}
}
\details{
A typical record in fast format is much like below:

>SPAC212.06c|SPAC212.06c.1|18558|18974|1||I|protein_coding| 
DNA helicase in rearranged telomeric region, truncated|
CTACACATTACGCTGAGAGGTAAAATACTCTGACAACATTCGTTCGATTGTATAAAACAA
AATCCAGCCGAAACGATTGTTGTCAGTAATCAAGATTACGATCTAAATTGAGTACCAAGA
CAAAACGAAATGGTTAAAAAGTTAAAGTCGTTTTTGTATGGACACAATTTCTATAAAATA
GACATGAGTAAAATCTCGCTATTTGTTTGTTATTGTGGAATAATGAAGAGTCATGGGAGA
TGAATGTTGTAAACGATGGCATAGAATTGGTAACGAAAAGTGAAATCGTTGGGATCAACT
ATTTCAGTATTTTGTTTAAAGAAAATGTTGAACTCGACAAGTAATGAGAGGTGGTGCTTT
CGTTAAATAATGAGTGGTGGTTACGGTTATACAGGATATGATATGTGTATGGTGAGA  
  
The header line includes items for  
    Systematic_ID
    UTR ID
	Gene Start  
	Gene_End  
	Strand  
	Chromosome  
	Feature Type  
	Description
}
\value{
A data frame with rows for UTR and columns for annotation items and sequence.
}
\references{
https://www.pombase.org/downloads/utr
}
\keyword{ methods }
