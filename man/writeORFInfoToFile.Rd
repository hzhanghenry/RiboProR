\name{writeORFInfoToFile}
\alias{writeORFInfoToFile}
\title{
Save ORF Information to File
}
\description{
Write ORF information to tab-delimited text file including column headers.
}
\usage{
writeORFInfoToFile(seq_info, start_codons, stop_codons, out_file)
}
\arguments{
  \item{seq_info}{
Data frame, UTR/CDS info including chromosome,  start and stop position,  
strand, sequence, and type (CDS, 5_UTR or 3_UTR).
}
  \item{start_codons}{
Character vector, start codons.
}
  \item{stop_codons}{
Character vector, stop codons.
}
  \item{out_file}{
Character vector, name (and path) of the file to write.
}
}
\details{
The output file will include, for each ORF, the gene_name, chromosome,   
seq_start, seq_stop, strand, orf_id, orf_start, orf_stop, start_codon,  
context, orf_length, dis_from_cap, dis_to_MUG.
}
\value{
None
}
\keyword{ methods }
