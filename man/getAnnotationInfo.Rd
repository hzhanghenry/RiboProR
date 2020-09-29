\name{getAnnotationInfo}
\alias{getAnnotationInfo}
\title{
Extract Annotation Information for a List of Genes
}
\description{
Extract gene ID, gene name, and description for a list of genes from annotation file.
}
\usage{
getAnnotationInfo(annotation, gene_list, id_column=1, name_column=2, description=3)
}
\arguments{
  \item{annotation}{
  A data frame with at least 3 columns for gene ID, gene names, and descriptions.
}
  \item{gene_list}{
  Character vector, list of gene names for which the annotation information is extracted.
}
  \item{id_column}{
  Positive integer, column of gene ID in annotation file.
}
  \item{name_column}{
  Positive integer, column of gene name in annotation file.
}
  \item{description}{
  Positive integer, column of gene description in annotation file.
}
}
\value{
A data frame with 3 columns (gene ID, gene name, and description) for the list of genes.
}
\keyword{ methods }

