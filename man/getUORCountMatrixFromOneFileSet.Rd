\name{getUORCountMatrixFromOneFileSet}
\alias{getUORCountMatrixFromOneFileSet}
\title{
Generate uORF Counts Matrix from One Set of Files
}
\description{
Generate uORF count matrix from files which contains both fp counts for uORF and cds regions.
}
\usage{
getUORCountMatrixFromOneFileSet(directory_name, 
	uorf_file_pattern, mrna_file_pattern=NULL, 
	count_column=2, rowname_column=0, has.header=TRUE, 
	uorf_ID_patter1 = "^.{4}-", uorf_ID_patter2=NULL)
}
\arguments{
  \item{directory_name}{
  Character vector, name of the directory where raw counts files are stored.
}
  \item{uorf_file_pattern}{
  Character vector, the common pattern of raw counts files.
}
  \item{mrna_file_pattern}{
  Character vector, the common pattern of RNASeq raw count files if using RNASeq reads counts instead of cds fp counts.
}
  \item{count_column}{
  Positive integer, the column of raw count in counts files.
}
  \item{rowname_column}{
  Positive integer, the column in raw counts file for row names of new matrix. set to 0 if raw counts files were saved by write.table() with row.names=TRUE.
}
  \item{has.header}{
  Logical, if the raw counts files have column headers.
}
  \item{uorf_ID_patter1}{
  Character vector, unique pattern used to extract gene ID from uorf ID.
}
  \item{uorf_ID_patter2}{
  Character vector, the second unique pattern, if exists, used to extract gene ID from uorf ID.
}
}
\value{
An integer matrix with raw footprints (fp) counts of both uORF and fp counts (or reads counts, if using RNASeq) cds regions. Rows of the matrix are for uORF/genes and columns for samples. The left half of the matrix is fp counts from cds regions of each gene and the right half is fp counts from uORF region of same gene. The samples in uORF part and cds part have same orders.
}
\keyword{ methods }

