#*
#*		File:	RiboProfiling_Default_Parameters.R
#*
#*		Functions to get default parameters for ribosomal footprinting 
#*		analsis. Default parameters are based on Dr. Nocholas Ingoliad's 
#*		source. 
#*
#*		Reference:
#*
#*		Nicholas T. Ingolia, Gloria A. Brar, Noam Stern-Ginossar, Michael S. 
#*		Harris, Gaelle J.S. Talhouarne, Sarah E. Jackson, Mark R. Wills,
#*		and Jonathan S. Weissman (2014). Ribosome Profiling Reveals Pervasive
#*		Translation Outside of Annotated Protein-Coding Genes. Cell Reports  
#*		8, 1365~1379.
#*
#*		Last revised on April 18, 2018
#*	________________________________________________________________________
#*	<parameters><parameters><parameters><parameters><parameters><parameters>




#	Default parameters
#	===========================================

getAllDefaultParameters <- function()
{
	parameters <- list(
		shift_start=100, shift_end=100,
		inset_5=34, inset_3=31,
		min_len=25, max_len=34,
		minLenFract=0.05,
		minStartRange=-17, maxStartRange=-8,
		minEndRange=-22, maxEndRange=-13,
		startShift=3,endShift=-2
	);

	return (parameters);
}

getDefaultInset <- function()
{
	return (list(utr5Inset3=6, cdsInset5=45, 
		cdsInset3=15, utr3Inset5=6));
}


getZeroInset <- function()
{
	return (list(utr5Inset3=0, cdsInset5=0, 
			cdsInset3=0, utr3Inset5=0));
}
	
getDefaultRiboSeqAsites <- function()
{
	return (data.frame(asite=c(14,14,14,14,15,15), 
				row.names=seq(26,31)));
}

getDefaultRNASeqAsites <- function()
{
	return (data.frame(asite=rep(15, 34), 
				row.names=seq(18,51)));
}
	
getDefaultFlank <- function()
{
	return (c(befor=-100, after=100));
}

getDefaultCdsBody <- function()
{	
	return (c(at_start=34, at_end=31));
}

getDefaultLengths <- function()
{
	return (c(min_len=25, max_len=34));
}
	
getDefaultMinLenFract <- function()
{
	return (minFraction=0.05);
}


getDefaultStartRange <- function()
{
	return (c(-17, -8))
}
getDefaultEndRange <- function()
{
	return (c(-22, -13))
}


getDefaultStartShift <- function()
{
	return (defaultStartShift=3);
} 
getDefaultEndShift <- function()
{
	return (defaultEndShift=-2);
}


getDefaultStartCodons <- function()
{
	start_codons <- c("ATG", "CTG", "GTG", "TTG", "AAG", 
					  "ACG", "AGG", "ATA", "ATC", "ATT"); 
	return (start_codons);
}
getDefaultStopCodons <- function()
{
	stop_codons <- c("TAA", "TAG", "TGA");
	return (stop_codons);
}
getDefaultExtraBounds <- function()
{
	return (c(100, 100));
}
getDefaultMinNumberOfAminoAcid <- function()
{
	return (2);
}


#*
#*	End of RiboProfiling_Default_Parameters.R
#*	________________________________________________________________________
#*	<parameters><parameters><parameters><parameters><parameters><parameters>
