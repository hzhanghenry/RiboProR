#*
#*    File: RiboProfiling_Footprint_Assignment.R
#*
#*    Calculate frame size for transcripts defined in annotation file
#*    (bed format) and for reads in bam file aligned to the transcripts
#*
#*    Source code from supplyment of reference:
#*
#*    Nicholas T. Ingolia, Gloria A. Brar, Noam Stern-Ginossar, Michael S.  
#*    Harris, Gaelle J.S. Talhouarne, Sarah E. Jackson, Mark R. Wills,
#*    and Jonathan S. Weissman (2014). Ribosome Profiling Reveals Pervasive
#*    Translation Outside of Annotated Protein-Coding Genes. Cell Reports  
#*    8, 1365~1379.
#*
#*    Function defined:
#*
#*    1.  alignASites()
#*    2.  getAlignASites()
#*    3.  getASiteProfile()
#*    4.  trxRegionCountSizes()
#*    5.  trxRegionCountAligns()
#*    6.  countAligns()
#*    7.  countSizes()
#*    8.  regionCountFrame()
#*
#*    Last revised on April 18, 2018
#*  ___________________________________________________________________________
#*  <fp assignment><fp assignment><fp assignment><fp assignment><fp assignment>


    
#*  =========================================================================
#*
#*    1.  alignASites()
#*
#*    Covert GAligments objects for a transcript to a GAlignments object
#*    having the A site nucleotides only.
#*
#*    Arguments:
#*
#*        asiteOffsets:  A data frame of one column with row names are 
#*                       fragment length and column is for A sites.
#*        alns:          A GAlignments objects for a transcript.
#*
#*    Return:  A GAlignments object for A site nucleotides only.
#*
#*
alignASites <- function(asiteOffsets, alns)  {
    alnLengths <- qwidth(alns)
    alnAOffs <- asiteOffsets[as.character(alnLengths), "asite"]

    alnDelta <- ifelse(strand(alns) == "-", alnLengths-alnAOffs, alnAOffs+1)
                       
    alnsGood <- alns[!is.na(alnDelta)]
    alnDelta <- alnDelta[!is.na(alnDelta)]

    if (length(alnsGood) > 0) {
        qnarrow(alnsGood, start = alnDelta, width = 1)
    } else {
        GRanges()
    }
}


#*	=========================================================================
#*
#*    2. getAlignASites()
#*
#*    Get aligned A sites for a transcript. This only call the function
#*    above (alignASites()).
#*
#*    Arguments:
#*
#*        asiteOffsets:  A data frame of one column with row names are 
#*                       fragment length and column is for A sites.
#*        bamfle:        character vector, bam file name and path.
#*        trx:           GRange object for a transcript.
#*
#*    Return:  A GAlignments object for A site nucleotides only.
#*
#*
getAlignASites <- function(asiteOffsets, bamfile, trx)  {
    alignASites(asiteOffsets, getAlignments(bamfile, trx))
}


#*  =========================================================================
#*
#*    3.  getASiteProfile()
#*
#*    Covert aligned A sites to a table for relative position of A sites
#*    in each transcript
#*
#*    Arguments:
#*
#*        asiteOffsets:	 A data frame of one column with row names are 
#*                       fragment length and column is for A sites.
#*        bamfle:        character vector, bam file name and path.
#*        trx:           GRange object for a transcript.
#*
#*    Return:		
#*
#*        Positive integer vector, a site status for each nucleotides
#*        of the transcript.
#*
#*
getASiteProfile <- function(asiteOffsets, bamfile, trx)  {
    asites <- getAlignASites(asiteOffsets, bamfile, trx)

    if (length(asites) > 0) { 
        tabulate( relativeWithin(start(asites), trx),
            nbins = sum(width(trx)) )
    } else {
        tabulate( as.integer(c()), nbins = sum(width(trx)) )
    }
}


#*  =========================================================================
#*
#*    4.  trxRegionCountSizes()
#*
#*    Calculate size of whole transcript, cds, 5'UTR, and 3'UTR of a 
#*    transcript in bed file.
#*
#*    Arguments:
#*
#*        insets:  List of integers, insets in nucleotides to avoid  
#*                 start and stop positions. 
#*        trx:     GRange object for transcripts, output from import()
#*        cds:     IRanges object for cds.
#*
#*    Return:   List of integers for sizes of cds, 5-UTR, and 3-UTR.
#*
#*
trxRegionCountSizes <- function(insets, trx, cds)  {

    sizes <- list( trx = sum(width(trx)) )
    
    if (length(cds) == 1)  {
        cdsQStart <- start(cds) + insets$cdsInset5
        cdsQEnd   <- end(cds) - insets$cdsInset3
        
        sizes$cds <- if (cdsQEnd < cdsQStart) { 
            NA 
        } else { 1 + cdsQEnd - cdsQStart }
        
        utr5QEnd   <- start(cds) - (insets$utr5Inset3 + 1)
        sizes$utr5 <- if (utr5QEnd < 0) { NA } else { 1 + utr5QEnd }

        utr3QStart <- end(cds) + (insets$utr3Inset5 + 1)
        sizes$utr3 <- if (utr3QStart >= sum(width(trx))) { 
            NA 
        } else { 1 + sum(width(trx)) - utr3QStart }
    } else {
        sizes$cds  <- NA
        sizes$utr5 <- NA
        sizes$utr3 <- NA
    }

    sizes
}


#* =========================================================================
#*
#*    5.  trxRegionCountAligns()
#*
#*    Check out if reads of a transcript has A site. If yes, get the total 
#*    counts for transcript, cds, 5-UTR, and 3-UTR. Otherwise, set counts 
#*    of the transcript to 0, and 0 or NA( if no cds defined) for cds, 5-UTR, 
#*    and 3-UTR.
#*
#*    Arguments:
#*		
#*        asiteOffsets:	 A data frame of one column with row names are 
#*                       fragment length and column is for A sites.
#*        insets:        List of integer, insets in nucleotides for    
#*                       avoiding start and stop. There are two build-in  
#*                       insets: default and zero insets.
#*        bamfile:       Character vector, bam file name and path.
#*        trx:           GRange object for a transcript.
#*        cds:           IRange object for cds in the transcript.
#*
#*    Return:  List of sizes for cds, 5-UTR, and 3-UTR of a transcript
#*
#*
trxRegionCountAligns <- function(asiteOffsets, insets, bamfile, trx, cds) {

    genomicAsites <- getAlignASites(asiteOffsets, bamfile, trx);

    if (length(genomicAsites) > 0) {

        asites <- relativeWithin(start(genomicAsites), trx);
        counts <- list( trx = sum(!is.na(asites)) );
        
        if ( length(cds) == 1 ) {
            cdsQStart <- start(cds) + insets$cdsInset5;
            cdsQEnd   <- end(cds) - insets$cdsInset3;
            
            counts$cds <- if (cdsQEnd < cdsQStart) { 
				NA 
            } else { 
                sum(asites >= cdsQStart & asites <= cdsQEnd, na.rm=TRUE) 
            }
            
            utr5QEnd <- start(cds) - (insets$utr5Inset3 + 1);
            counts$utr5 <- if (utr5QEnd < 0) { 
                NA 
            } else { 
                sum(asites <= utr5QEnd, na.rm=TRUE) 
            }
            
            utr3QStart  <- end(cds) + (insets$utr3Inset5 + 1);
            counts$utr3 <- if (utr3QStart >= sum(width(trx))) { 
                NA 
            } else { 
                sum(asites >= utr3QStart, na.rm=TRUE) 
            }
        } else {
            counts$cds  <- NA
            counts$utr5 <- NA
            counts$utr3 <- NA           
        }
    } else {
        counts <- list( trx = 0 )

        if ( length(cds) == 1) {
            counts$cds <- 0
            counts$utr5 <- 0
            counts$utr3 <- 0
        } else {
            counts$cds <- NA
            counts$utr5 <- NA
            counts$utr3 <- NA
		}
    }
    
    counts
}


#*  =========================================================================
#*
#*    6.  countAligns()
#*
#*    Count the reads in a bam file which hits transcript, cds, 5-UTR,
#*    and 3-UTR for all transcripts. 
#*
#*    Arguments:
#*		
#*        asiteOffsets:  A data frame of one column with row names are 
#*                       fragment length and column is for A sites.
#*        insets:        List of integer, insets in nucleotides for    
#*                       avoiding start and stop position. There are   
#*                       two build-in insets: default and zero insets.
#*        bamfile:       Character vector, bam file name and path.
#*        trxBed:        GRange object for all transcripts.
#*
#*    Return:	
#*
#*        List of positive integers, counts for whole transcript, cds,  
#*        5-UTR, and 3-UTR of all transcripts.  
#*
#*
countAligns <- function(asiteOffsets, insets, bamfile, trxBed) {

    ttl <- length(trxBed)
    countOne <- function(i) {
        if (i %% 100 == 0) { 
            print(c(i, ttl, as.character(trxBed$name[i]))) 
        }
        trxGR <- transcriptGRanges(trxBed[i])
        cdsIR <- transcriptCdsIRanges(trxGR, 
            start(trxBed$thick)[[i]], 
            end(trxBed$thick)[[i]]);
        trxRegionCountAligns(asiteOffsets, insets, bamfile, trxGR, cdsIR);
    }
    counts <- mclapply(seq(1,ttl), countOne)
    names(counts) <- trxBed$name
    counts
}


#*  =========================================================================
#*
#*    7.  countSizes()
#*
#*    Count transcript sizes for each transcript in bed file. Output will
#*    be used for GetRegionCountFrame()
#*		
#*    Arguments:
#*
#*        insets:  list of integer, insets in nucleotides for avoiding   
#*                 star tand stop.    
#*        trxBed:  GRange object of all transcripts.
#*
#*    Return: 
#*
#*        List of sizes of cds, 5-UTR, and 3-UTR of all transcript.
#*
#*
countSizes <- function(insets, trxBed) {
    countOne <- function(i) {
        trxGR <- transcriptGRanges(trxBed[i])
        cdsIR <- transcriptCdsIRanges(trxGR, start(trxBed$thick)[[i]], 
                    end(trxBed$thick)[[i]]);
        trxRegionCountSizes(insets, trxGR, cdsIR);
    }    
    sizes <- mclapply(seq(1,length(trxBed)), countOne);
    names(sizes) <- trxBed$name;
    sizes
}


#*  =========================================================================
#*
#*    8.  regionCountFrame()
#*
#*    Covert lists for sizes of transcript, cds, 5-UTR, and 3-UTR to
#*    a data frame.
#*
#*    Argument:	
#*
#*        List of integer lists, read counts in cds, 5-UTR, and 3-UTR of   
#*        each transcript.
#*
#*    Return: 
#*
#*        A data frame with rows for transcript names and columns for size
#*        of transcript, cds, 5-UTR, and 3-UTR.
#*
#*
regionCountFrame <- function(counts)  {
    data.frame(trx  = unlist(lapply(counts, function(r) { r$trx })),
               cds  = unlist(lapply(counts, function(r) { r$cds })),
               utr5 = unlist(lapply(counts, function(r) { r$utr5 })),
               utr3 = unlist(lapply(counts, function(r) { r$utr3 })),
               row.names = names(counts));
}


#*
#*  End of RiboProfiling_Footprints_Assignment.R
#*  ___________________________________________________________________________
#*  <fp assignment><fp assignment><fp assignment><fp assignment><fp assignment>
