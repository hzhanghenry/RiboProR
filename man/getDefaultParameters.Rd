\name{getDefaultParameters}
\alias{getAllDefaultParameters}
\alias{getDefaultInset}
\alias{getZeroInset}
\alias{getDefaultRiboSeqAsites}
\alias{getDefaultRNASeqAsites}
\alias{getDefaultFlank}
\alias{getDefaultCdsBody}
\alias{getDefaultLengths}
\alias{getDefaultMinLenFract}
\alias{getDefaultStartRange}
\alias{getDefaultEndRange}
\alias{getDefaultStartShift}
\alias{getDefaultEndShift}
\alias{getDefaultStartCodons}
\alias{getDefaultStopCodons}
\alias{getDefaultExtraBounds}
\alias{getDefaultMinNumberOfAminoAcid}
\title{
Get Default Parameters
}
\description{
Methods used for get default parameters.
}
\usage{
getAllDefaultParameters()
getDefaultInset()
getZeroInset()
getDefaultRiboSeqAsites()
getDefaultRNASeqAsites()
getDefaultFlank()
getDefaultCdsBody()
getDefaultLengths()
getDefaultMinLenFract()
getDefaultStartRange()
getDefaultEndRange()
getDefaultStartShift()
getDefaultEndShift()
getDefaultStartCodons()
getDefaultStopCodons()
getDefaultExtraBounds()
getDefaultMinNumberOfAminoAcid()
}
\value{
getAllDefaultParameters() return a list of integers used for ribosomal footprints framing:
 		
	shift_start:100 
	shift_end:100
	inset_5: 34 
	inset_3: 31
	min_len: 25 
	max_len: 34
	minLenFract: 0.05,
	minStartRange: -17, 
	maxStartRange: -8,
	minEndRange: -22, 
	maxEndRange: -13,
	startShift: 3,
	endShift: -2
  
getDefaultInset() return a list of integers in order to avoid start and stop codons:

	utr5Inset3: 6, 
	cdsInset5: 45, 
	cdsInset3: 15, 
	utr3Inset5: 6
		
getZeroInset() return a list of 4 zeros for utr5Inset3, cdsInset5, cdsInset3, and utr3Inset5.  

getDefaultRiboSeqAsites() return a data frame of 1 colum for a-sites for Riboseq reads with length 26:31.  

getDefaultRNASeqAsites() return a data frame of 1 colum for a sites for RANSeq reads with length 18:51.  

getDefaultFlank() return an integer vector of length 2 for flank regions at both end of gene cds: -100 and 100.  

getDefaultCdsBody() return an integer vector of length 2 to calculate CDS body range: 34 and 31.  

getDefaultLengths() return an integer vector of length 2 for minimum and mazimum read length: 25 and 34.  

getDefaultMinLenFract() return a numeric for minimum required fraction of read length: 0.05.  

getDefaultStartRange() return an integer vector of length 2 for range before start codon.  

getDefaultEndRange() return an integer vector of length 2 for range before stop codon.  

getDefaultStartShift() return a integer (3) for start position when shift starts.  

getDefaultEndShift() return a integer (-2) for stop position when shif stops.  

getDefaultStartCodons() return a character vector of length 10 for start codons.  

getDefaultStopCodons() return a character vector of length 3 for stop codons.  

getDefaultExtraBounds() return an integer vector of length 2 for extra bound of a region (such as cds).  

getDefaultMinNumberOfAminoAcid() return an integer (2) for minimum number of ammino acids required for a gene.  
}
\keyword{ methods }

