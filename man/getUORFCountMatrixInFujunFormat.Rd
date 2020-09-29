\name{getUORFCountMatrixInFujunFormat}
\alias{getUORFCountMatrixInFujunFormat}
\title{
Generate uORF FP Counts Matrix in Fujun Format
}
\description{
Generate a uORF footprints counts matrix in a special format. See Details.
}
\usage{
getUORFCountMatrixInFujunFormat(file_directory, 
	uorf_file_pattern = "uorf_conserved_fpcounts.txt", 
	cds_file_pattern = "cds_for_uorf_fpcounts.txt", 
	count_column=2, rowname_column=0, has.header=TRUE)
}
\arguments{
  \item{file_directory}{
  Character vector, name of a directory where raw footprints counts files are stored.
}
  \item{uorf_file_pattern}{
  Character vectors, common pattern in uORF footprints counts files.
}
  \item{cds_file_pattern}{
  Character vectors, common pattern in cds footprints counts files.
}
  \item{count_column}{
  Possitive interger, the column of raw counts in counts files.	
}
  \item{rowname_column}{
  Possitive interger, which column will be used as row names. Set to 0 if count files are output from R write.table() with row.names=TRUE. 
}
  \item{has.header}{
  Logic, if the input files have headers.
}
}
\details{
The Fujun Format is specifically designed for relative risobomal occupancy (RRO, uORF fp counts divided by cds fo counts) analysis. Since in most cases the uORF fp counts is much lower than the counts of cds region in same gene, simply normalizing the uORF counts with relative cds fp counts together may extremely shift the values of log2 change ratio. With Fujun Format, the coumns for uORF counts in the matrix contains fp counts for all uORFs and all cds regions, the columns for cds region fp counts have cds region fp counts of genes matched to each uORF plus fp counts of all genes. A default DESeq2 analysis could be applied with data matrix in this format.
}
\value{
An integer matrix with raw footprints (fp) counts of both uORF and cds regions. Rows of the matrix are for uORF/genes and columns for samples. The left half of the matrix is fp counts from uORF regions of each gene and the right half is
fp counts from cds region of same gene. The samples in uORF part and cds part have same orders.
}
\note{
The input files are outputs based on two different BED files but the rows in two BED files are matched each other.
The uORF ID/names in raw uORF counts file must contain gene name and uORF ID separated by the 'seperator', e.g., "YGL134W.255638" where YGL134W is gene name, "." is separator, and 255638 is uORF ID.
}
\keyword{ methods }
