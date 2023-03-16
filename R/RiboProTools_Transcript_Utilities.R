#*    File: RiboProfiling_Transcript_Utilities.R
#*
#*    Functions to manipulate transcripts with GenomicRanges package.
#*    Most source code are from Ingolia's Lab and were Modified for 
#*    easy understanding.
#*
#*    Reference:
#*
#*    Nicholas T. Ingolia, Gloria A. Brar, Noam Stern-Ginossar, Michael 
#*    S. Harris, Gaelle J.S. Talhouarne, Sarah E. Jackson, Mark R. Wills,
#*    and Jonathan S. Weissman (2014). Ribosome Profiling Reveals Pervasive
#*    Translation Outside of Annotated Protein-Coding Genes. Cell Reports
#*    8, 1365-1379.
#*
#*    Function defined:
#*
#*    1.   relativeWithin()
#*    2.   absoluteOuter()
#*    3.   irangeOuter()
#*    4.   transcriptGRanges()
#*    5.   transcriptCdsIRanges()
#*    6.   bedCdsIRangesList()
#*    7.   bedGRangesList()
#*    8.   trxRegions()
#*    9.   getAllAlignments()
#*    10.  getAlignments()
#*    11.  getSequence()
#*    12.  getBamFlagStat()
#*    13.  ReadBamToTable()
#*
#*    Last revised on April 18, 2018
#*  ________________________________________________________________________
#*  <transcript-utils><transcript-utils><transcript-utils><transcript-utils>




#*  =========================================================================
#*
#*    1.  relativeWithin()
#*
#*    Convert an absolute genomic query position (qpos) to a transcript-
#*    relative coordinate position within a transcript whose coordinates 
#*    are given by a GRanges of exons (outer).  NOTE that the exons in
#*    outer are in ascending numeric order regardless of the strand.
#*		
#*    Arguments:
#*
#*        qpos:   A positive integer, a genomic coordinate inside of a 
#*                transcript.							
#*        outer:  GRanges object, genomic positions of a transctipt.
#*
#*    Return:
#*
#*        Positive integer, a poistion relative to the start of 
#*        transctipt or NA if the absolute position is not in any 
#*        exon (no hit).
#*
relativeWithin <- function(qpos, outer)  {

    outerOffsets <- start(outer);
    withinOffsets <- head(c(0, cumsum(width(outer))), -1);

    startHits <- findOverlaps(IRanges::IRanges(start=qpos, width=1), 
        ranges(outer), select="first");
  
    within <- (qpos - outerOffsets[startHits]) + withinOffsets[startHits];

    if (as.vector(strand(outer))[1] == "-") {
        within <- (sum(width(outer)) - 1) - within;
    }
    within
}


#* =========================================================================
#*
#*    2.  absoluteOuter()
#*
#*    Convert a transcript-relative query position (qpos) to an absolute
#*    genomic coordinate based on a transcript whose coordinates are given
#*    by a GRanges of exons (outer). NOTE that the exons are in ascending 
#*    numeric order regardless of the strand.
#*
#*    Arguments:
#*
#*        qpos:   Positive integer, a position relative to the start 
#*                of a transcript.
#*        outer:  GRanges object, genomic positions of a transctipt.
#*
#*    Return:
#*
#*        Positive integer, geneomic coordinate of the transcript-relative 
#*        position.
#*
#*
absoluteOuter <- function(qpos, outer)  {

    outerOffsets <- start(outer);
    withinOffsets <- head(c(0, cumsum(width(outer))), -1);
    inner <- IRanges::IRanges(start=withinOffsets, width=width(outer));

    if (as.vector(strand(outer))[1] == "-") {
        qpos <- (sum(width(outer))-1) - qpos;
    }

    startHits <- findOverlaps(IRanges::IRanges(start=qpos, width=1), 
        inner, select="first");

    (qpos - withinOffsets[startHits]) + outerOffsets[startHits];
}


#*  =========================================================================
#*
#*    3.  irangeOuter()
#*
#*    Convert transcript-relative query IRanges (qranges) to an absolute
#*    genomic GRanges based on a transcript whose coordinate are given by
#*    a GRanges of exons (outer).
#*
#*    Arguments:
#*
#*        qranges:  IRanges object, a relative positions to a start 
#*                  of transcript in GRanges format.
#*        outer:    GRanges object, genomic positions of a transctipt.
#*
#*    Return:  a Granges object.
#*
#*    Note: The relative range is in the second exon of transcript.
#*

irangeOuter <- function(qranges, outer) {

    starts <- absoluteOuter(start(qranges), outer);
    ends <- absoluteOuter(end(qranges), outer);

    if (as.vector(strand(outer))[1] == "-") {
        masks <- GRanges(seqnames = as.vector(seqnames(outer))[[1]],
                    ranges = IRanges::IRanges( start=ends, end=starts ),
                    strand = "-" );
    } else {
        masks <- GRanges(seqnames = as.vector(seqnames(outer))[[1]],
                  ranges = IRanges::IRanges( start=starts, end=ends ),
                  strand = "+" );  
    }

    intersect(masks, outer);
}


#* =========================================================================
#*
#*    4.  transcriptGRanges()
#*
#*    bedGRange should be a 1-element GRanges slice. The metadata columns
#*    should include a "blocks" entry as from an rtracklayer BED file.
#*
#*    A GRanges object is returned with one range per exon.
#*    NOTE that the exons are in ascending numeric order regardless of the 
#*    transcript strand.
#*
#*    Argument:
#*
#*        bedGRange:  a GRange object from import() in rtracklayer package.
#*
#*    Return: GRanges object with one range per exon.
#*
#*
transcriptGRanges <- function(bedGRange)  {

    blocks <- mcols(bedGRange)$blocks[[1]];
    exons <- GenomicRanges::shift(blocks, shift = start(bedGRange)[[1]] - 1);
	
    sn <- seqnames(bedGRange)[1];
    #runLength(sn) <- length(exons);
	
    st <- strand(bedGRange)[1];
    #runLength(st) <- length(exons);
    
    GRanges(seqnames = sn,
            ranges = exons,
            strand = st);
}


#*  =========================================================================
#*
#*    5.  transcriptCdsIRanges()
#*
#*    trx is a transcript GRanges wiht one range per exon, and thickStart
#*    & thickEnd are the genomic start and end of the coding sequence,
#*    e.g. taken from start(...) and end(...) on the $thick IRanges from a
#*    bed file.
#*
#*    Arguments:
#*
#*        trx:         GRanges object for a transcript.
#*        thickStart:  positive integer, start position of a CD.
#*        thickEnd:    positive integer, end position of a CD.
#*
#*    Return:  An IRanges object for start and end position of a CD.
#*
transcriptCdsIRanges <- function(trx, thickStart, thickEnd) {

    startWithin <- relativeWithin(thickStart, trx);
    endWithin <- relativeWithin(thickEnd, trx);

    if (is.na(startWithin) || is.na(endWithin)) {
        IRanges::IRanges();
    } else {
        IRanges::IRanges(start=min( startWithin, endWithin ), 
                end=max( startWithin, endWithin ) );
    }
}


#* =========================================================================
#*
#*    6.  bedCdsIRangesList()
#*
#*    For each transcript in a BED file and the associated transcript
#*    GRanges as computed by bedGRangesList, find the transcript-local 
#*    CDS as per transcriptCdsIRanges and collect these as an IRangesList.
#*
#*    Arguments:
#*
#*        bed:     GRanges object, such as output from import().
#*        bedgrl:  GRanges list, output of bedGRangesList().
#*
#*    Return:  IRanges list for each exon (CDS).
#*
bedCdsIRangesList <- function(bed, bedgrl)  {
    irl <- IRanges::IRangesList(mclapply(1:length(bed), 
                function(i) { 
                    transcriptCdsIRanges(bedgrl[[i]], 
                    start(bed$thick)[[i]], 
                    end(bed$thick)[[i]] ) } ));

    names(irl) <- bed$name;
    irl;
}


#* =========================================================================
#*
#*    7. bedGRangesList()
#*
#*    For each transcript in a BED file GRanges, get the GRanges of its
#*    exons as per transcriptGRanges and collec these into a GRangesList
#*    with names taken from the name metadata column of the BED file
#*    GRanges.
#*
#*    Argument:
#*
#*        bed:  GRangesobject, such as the output of import().
#*
#*    Return:  A GRanges list object.
#*
bedGRangesList <- function(bed) {
    grl <- GRangesList(
                mclapply(1:length(bed), 
                function(i) { transcriptGRanges(bed[i]) })
            );
    names(grl) <- bed$name;
    grl;
}


#* =========================================================================
#*
#*    8.  trxRegions()
#*
#*    Calculate a GRangesList with named GRanges giving genomic  
#*    coordinates of transcript regions: trx for the entire transcript,  
#*    cds for the coding sequence, utr5 for the 5' UTR, and utr3 for 
#*    the 3' UTR.
#*
#*    Arguments:
#*
#*        trx: 	GRanges object for a transcript.
#*        cds:	GRanges object for cds in a transcript.
#*
#*        insetUtr5Start, insetUtr5End,	
#*        insetCdsStar, insetCdsEnd,
#*        insetUtr3Start, insetUtr3End: 
#*
#*            positive integer, for adjustment of cds, 5-UTR, and 3-UTR 
#*            positions.
#*
#*    Returned:	
#*
#*        GRanges list represent cds, 5-UTR, and 3-UTR for a transcript.
#*
#*
trxRegions <- function(trx, cds,
    insetUtr5Start=0, insetUtr5End=0,
    insetCdsStart=0,  insetCdsEnd=0,
    insetUtr3Start=0, insetUtr3End=0)  {

    rgns <- list( trx = trx );

    if (length(cds) == 1) {
        utr5Start <- insetUtr5Start;
        utr5End   <- (start(cds) - 1) - insetUtr5End;
        cdsStart  <- start(cds) + insetCdsStart;
        cdsEnd    <- end(cds) - insetCdsEnd;
        utr3Start <- (end(cds) + 1) + insetUtr3Start;
        utr3End   <- (sum(width(trx)) - 1) - insetUtr3End;

        if (cdsStart < cdsEnd) {
            rgns$cds <- irangeOuter(IRanges::IRanges(start=cdsStart, 
                end=cdsEnd), trx);
        }

        if (utr5Start < utr5End) {
            rgns$utr5 <- irangeOuter(IRanges::IRanges(start=utr5Start, 
                end=utr5End ), trx);
        }

        if (utr3Start < utr3End) {
            rgns$utr3 <- irangeOuter(IRanges::IRanges(start=utr3Start, 
                end=utr3End), trx );
        }
    }

    rgns
}


#*  =========================================================================
#*
#*    9.  getAllAlignments()
#*
#*    GAlignments containing all reads overlapping the genomic extent 
#*    of a primary transcript on the correct strand, including unspliced 
#*    and purely intronic reads.
#*
#*    Note: readGAlignmentsFromBam has been deprecated. Instead, use the
#*          function readGAlignments().
#*
#*    Arguments:
#*
#*        bamfile:  Character vector, the name (and path) of the bam 
#*                  file to be read.
#*        trx:      GRanges object for a transcript.
#*
#*    Return:
#* 
#*        A GAlignments object with all reads overlapping with 
#*        the transcript.
#*
#*
getAllAlignments <- function(bamfile, trx)  {

    range <- GRanges( seqnames = seqnames(trx)[1],
            ranges = IRanges::IRanges(start = min(start(trx)), 
						end = max(end(trx)))
    );
    
    isminus <- (as.vector(strand(trx))[1] == "-")
    
    readGAlignments(bamfile, param=ScanBamParam(which=range, 
        flag=scanBamFlag(isMinusStrand=isminus)))
}


#* =========================================================================
#*
#*    10. getAlignments()
#*
#*    GAlignments containing all reads compatible with a processed 
#*    transcript on the correct strand, as per findSpliceOverlaps.
#*
#*    Arguments:
#*
#*    bamfile:  character vector, bam file name and path.
#*    trx:      GRanges object for a transcript.
#*
#*    Return:  A GAlignments object.
#*
getAlignments <- function(bamfile, trx) {

    alns <- getAllAlignments(bamfile, trx);
    
    spliceHits <- findSpliceOverlaps(alns, GRangesList(trx));
    compatibles <- spliceHits[mcols(spliceHits)$compatible];
    
    alns[queryHits(compatibles)]
}


#* =========================================================================
#*
#*    11.  getSequence()
#*
#*    Get sequence from fasta file for a transcript. 
#*
#*    Arguments:
#*
#*        fafile:  Character vector, the fasta file name (and path).
#*        trx:     GRanges object for a transcript.
#*
#*    Return: Character vector, DNA sequence of the transcript.
#*
getSequence <- function(fafile, trx)  {
    gseq <- unlist(scanFa(fafile, trx));

    if (as.vector(strand(trx))[1] == "-") {
        reverseComplement(gseq)
    } else {
        gseq
    }
}


#*  =========================================================================
#*
#*    12.  getBamFlagStat()
#*
#*    Scan bam flags to get statistics of the reads. This need samtools
#*    available from system (either path to samtools is included in
#*    user's PATH variable or the module has been loaded in HPC system. 
#*
#*    Arguments:
#*
#*        bamFiles:  Character vector, names of bam files.
#*        outFile:  Character vector, names of output file.
#*
#*    Return:  None, result is writen to file.
#*
getBamFlagStat <- function(bamFiles, outFile) {

    counts <- mclapply(bamFiles, function(bamFile) {
        res <- system(sprintf("samtools flagstat %s", bamFile), intern=TRUE)
        c(bamFile, res)
    })
    
    cat(do.call(c, counts), file=outFile, sep="\n")
}


#*  ===========================================================================
#*
#*    13.  Read a bam file and convert the object (list of list) to a 
#*         data frame
#*
#*    Reference:
#*
#*    Morgan M, Pages H, Obenchain V and Hayden N (2016). Rsamtools: 
#*    Binary alignment (BAM), FASTA, variant call (BCF), and tabix 
#*    file import. R package version 1.26.1, http://bioconductor.org/
#*    packages/release/bioc/html/Rsamtools.html.
#*
#*    Argument:  
#*
#*        bam.file: character vector, bam file name.
#*
#*    Return:    data frame with all contents from bam file.
#*
#* 
ReadBamToTable <- function(bam.file) {

    bam <- scanBam(bam.file);
    bam.field <- names(bam[[1]]);
  
    #*  a function to collapsing the list of lists into a single list
    #* 	-------------------------------------------------------------
    .unlist <- function (x)
    {
        x1 <- x[[1L]];
        if (is.factor(x1)){
            structure(unlist(x), class="factor", levels=levels(x1))
        } else {
            do.call(c, x)
        }
    }

    list <- lapply(bam.field, function(y) .unlist(lapply(bam, "[[", y)))
    bam.df <- do.call("data.frame", list)
    names(bam.df) <- bam.field

    #*  return a list that can be called as a data frame
    #*  ------------------------------------------------
    return(bam.df);
}


#*  =================================================================
#*
#*    14.  getTranscriptLength()
#*
#*    Function to get transcript length. The main purpose of this
#*    function is to get transcripts length with intron removed.
#*
#*    Argument:
#*
#*        transcripts: A GRange list with all transcripts
#*
#*    Return:
#*
#*        A data frame of 1 columns for transcript length and with 
#*        gene names as row names.
#*
getTranscriptLength <- function(transcripts) {
    counts <- mclapply(1:length(transcripts), 
        function(i) {sum(width(transcripts[i]$blocks[[1]]))	}
    )
    counts <- unlist(counts);
    return (data.frame(counts, row.names=transcripts$name))
}


#*
#*  End of RiboProfiling_Transcript_Utilities.R
#*  ________________________________________________________________________
#*  <transcript-utils><transcript-utils><transcript-utils><transcript-utils>


