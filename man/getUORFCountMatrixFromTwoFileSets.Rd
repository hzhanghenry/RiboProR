\name{getUORFCountMatrixFromTwoFileSets}
\alias{getUORFCountMatrixFromTwoFileSets}
\title{
Generate uORF Counts Matrix from Two Counts Files
}
\description{
Generate raw uORF counts matrix from two sets of counts files: 
footprint counts for uORF regions and footprints counts for 
gene cds regions.  
}
\usage{
getUORFCountMatrixFromTwoFileSets(file_directory, 
	uorf_file_pattern="uorf_conserved_fpcounts.txt", 
	cds_file_pattern="cds_for_uorf_fpcounts.txt", 
	count_column=2, rowname_column=0, has.header=TRUE, 
	seperator = "\\\.")
}
\arguments{
  \item{file_directory}{
Character vector, name of a directory where raw counts files 
are stored.  
}
  \item{uorf_file_pattern}{
Character vectors, common pattern in uORF counts files.  
}
  \item{cds_file_pattern}{
Character vectors, common pattern in cds counts files.  
}
  \item{count_column}{
Possitive interger, the column of raw counts in counts file.  	
}
  \item{rowname_column}{
Possitive interger, which column will be used as row names. Set 
to 0 if count files are output from R write.table() with 
row.names=TRUE.  
}
  \item{has.header}{
Logic, if the input files have headers.  
}
  \item{seperator}{
Character vector, unique pattern to separate gene name and uorf 
ID in uorf names. It CANNOT be any character used in gene name 
or uorf id.  
}
}
\value{
An integer matrix with raw footprints (fp) counts of both uORF 
and cds regions. Rows of the matrix are for uORF/genes and 
columns for samples. The left half of the matrix is fp counts 
from cds regions of each gene and the right half is fp counts 
from uORF region of same gene. The samples in uORF part and 
cds part have same orders.  
}
\note{
The two file sets have different rows, one is for uORF fp counts 
only and another is for cds region fp counts only.  
The uORF ID/names in raw uORF counts file must contain gene name 
and uORF ID separated by the 'seperator', e.g., "YGL134W.255638" 
where YGL134W is gene name, "." is separator, and 255638 is uORF 
ID.  
The data matrix generated with this function is for DESeq2 
analysis with normalization by subset, i.e., normalizing uORF 
and cds fp counts seperately. It MAY NOT be suitable for default 
normalization.  
}
\keyword{ methods }
