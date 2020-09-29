\name{filterCountMatrix}
\alias{filterCountMatrix}
\title{
Filter Data Matrix Based on RiboSeq fp Counts and/or RNAseq Counts
}
\description{
Filter data matrix based on RiboSeq fp counts and RNAseq counts. 
Any rows with mean of RiboSeq fp counts or mean of RNASeq counts 
below the definded values will be removed.  
}
\usage{
filterCountMatrix(count_matrix=NULL, mRNA_col=NULL, 
    ribo_col=NULL, mRNA_level=10, ribo_level=0)
}
\arguments{
  \item{count_matrix}{
Data matrix with both RiboSeq fp counts and RNASeq reads counts.  
}
  \item{mRNA_col}{
Positive integer, columns in the matrix for RNASeq reads counts.  
}
  \item{ribo_col}{
Positive integer, columns in the matrix for RiboSeq reads counts.  
}
  \item{mRNA_level}{
Numeric, threshold to filter rows by mean of RNASeq reads counts.  
}
  \item{ribo_level}{
Numeric, threshold to filter rows by mean of RiboSeq fp counts.  
}
}
\value{
A data matrix same as input matrix with unqualified rows removed.  
}
\examples{
data(ribo_pro_data);

mRNA_col <- c(4:6, 10:12); 
ribo_col <- c(1:3, 7:9);

mRNA_level <- 10;
Ribo_level <- 0;

count_matrix <- filterCountMatrix(ribo_pro_data, mRNA_col, 
					ribo_col, mRNA_level, Ribo_level);
}
\keyword{ methods }

