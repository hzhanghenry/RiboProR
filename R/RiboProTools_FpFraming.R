#*
#*    File: FpFraming.R
#*
#*    Base on the Haskell source code from Dr. Nocholas Ingoliad' lab 
#*    (https://github.com/ingolia-lab/RiboSeq)
#*
#*    Functions implemented:
#*
#*    1.   FpFramming()
#*    2.   getAnnotationFromBedFile()
#*    3.   getAlignmentsFromBamFile()
#*    4.   countingMetagenePosition()
#*    5.   countingFrames()
#*    6.   getPositionProfileTable()
#*    7.   initializeMetageneTable()
#*    8.   getFrameTable()
#*    9.   getASiteTable()
#*    10.  getReadPeaks()
#*    11.  getBestFrame()
#*    12.  validateParameters()
#*    13.  getMetageneFrames()
#*    14.  summerizeMetageneFrames()
#*    15.  initializeMetageneFrameTable()
#*
#*    Revised on April 18, 2018
#*    Last revised on August 6, 2018
#*  ________________________________________________________________________
#*  <fp Framing><fp Framing><fp Framing><fp Framing><fp Framing><fp Framing>




#*  ========================================================================
#*
#*    1.  FpFraming()
#*
#*    The main function to perform framing analysis based on Ingolia's 
#*    Haskell modules. Six tables will be generated and save to files
#*    if the save.file argument is set to TRUE. These tables are:
#*
#*    Two tables contain a 2-D metagene analysis of footprints  
#*    around the start and the end of protein-coding genes.
#*
#*    Two tables contain reading frame, totoal counts of reads of 
#*    each metagene position (relative to the start and the end of 
#*    protein-coding region). These two table could be used for
#*    metagene frame plot.
#*
#*    One table contains total numbers and fractions of all three  
#*    reading frames for each fragment length. This file could be 
#*    used to plot frame distributions (barplot). 
#*
#*    one table contains a table of a-sites for all fragment length.
#*    This table will be returned for downstream analysis.
#*
#*    Arguments:
#*
#*        bam_file:    Character vector, bam file name (and path).			
#*        bed_file:    Character vector, bed file name (and path).
#*        parameters:  List of integers for FpFraming parameters.
#*        save.file:   Logic, if all tables will be saved to files.
#*
#*    Return:
#*
#*        A data frame containing a-site for fragments with selected length. 
#*
FpFraming <- function(bam_file, bed_file, parameters=NULL, 
    validate.parameters=TRUE, save.file=TRUE) {

    #*  Arguments validation
    #*  ----------------------------------
    if(is.null(parameters)) {
        parameters <- getAllDefaultParameters()
    } else { 
        if(validate.parameters == TRUE)
            validateParameters(parameters);
    }

    shift_start <- parameters$shift_start; 
    shift_end	<- parameters$shift_end;
    inset_5 <- parameters$inset_5; 
    inset_3 <- parameters$inset_3;
    min_len <- parameters$min_len; 
    max_len <- parameters$max_len;

    #*  Gene annotation from bed file
    #*  ---------------------------------------------
    message("Reading ", bed_file, "...", appendLF = FALSE);	
    annotation <- getAnnotationFromBedFile(bed_file);
    cds_ranges  <- annotation[[2]];
    message("done!");

    #*  Scan all reads in bam file and keep reads in range of
    #*  min_len:max_len only
    #*  -----------------------------------------------------		
    message("Reading ", bam_file, "...", appendLF = FALSE);
    alignments <- getAlignmentsFromBamFile(bam_file);	
    alignments <- alignments[width(alignments) %in% min_len:max_len];
    message("done!");	

	
    #*  Get distance to start and end position of genes
    #*  -----------------------------------------------------	
    message("Counting metagene ...", appendLF = FALSE);	
    posProfile <- countingMetagenePosition(alignments, annotation,
        shift_start, shift_end);
    message("done!");

    message(paste("Getting position tables ... "), appendLF = FALSE);	
    at_start <- getPositionProfileTable(posProfile, at_which=1, 
        shift_start*-1, shift_end);
    at_end <- getPositionProfileTable(posProfile, at_which=2, 
        shift_start*-1, shift_end);			
    message("done!");

    #*  Get frame count for each reads
    #*  -----------------------------------------------------
    message("Counting frames ... ", appendLF = FALSE);		
    frameProfile <- countingFrames(alignments, cds_ranges, inset_5, inset_3);
    frame_table <- getFrameTable(frameProfile);
    message("done!");	
		
    #*  A-site table
    #*  ----------------------------------------------------
    message("Counting a-site ... ", appendLF = FALSE);	
    asite_table <- getASiteTable(frame_table, at_start, at_end, parameters);
    message("done.")
	
    if(save.file == TRUE){
	
        in_file_name <- basename(bam_file);
		
        message("Save at_start profile.");
        start_out_file <- sub("bam", "start_pos_len.txt", in_file_name)
        write.table(at_start, file=start_out_file, sep="\t",
            quote=FALSE, row.names=TRUE, col.names=TRUE);
	
        message("Save at_sop profile.");
        end_out_file <- sub("bam", "end_pos_len.txt", in_file_name);
        write.table(at_end, file=end_out_file, sep="\t",
            quote=FALSE, row.names=TRUE, col.names=TRUE);			

        message("Save frame profile.");
        frame_file <- sub("bam", "frame_len.txt", in_file_name);
        write.table(frame_table, file=frame_file, sep="\t",
            quote=FALSE, row.names=TRUE, col.names=TRUE);	

        message("Save asite profile.");			
        asite_file <- sub("bam", "asite_report.txt", in_file_name);
        write.table(asite_table, file=asite_file, sep="\t",
            quote=FALSE, row.names=TRUE, col.names=TRUE);	

        message("Save metagene frame profile");
        getMetageneFrames(posProfile, bam_file);
    }
	
    rm(alignments);
    rm(posProfile);
    rm(at_start);
    rm(at_end);
    rm(frameProfile);
    rm(frame_table);
    rm(annotation);
    rm(cds_ranges);
	
    return (asite_table);
}


#*  ===========================================================================
#*
#*    2.  getAnnotationFromBedFile()
#*
#*    Load BED file and get cds start and end positions for all genes
#*
#*    Argument:
#*
#*        bedFile: Character vector, file name (and path) for 
#*                 annotations in BED format. The file must be 
#*                 in Bed12 format (with thick columns defined).
#*
#*    Return:	
#*	
#*        List of two GRanges objects.  The first one contains all   
#*        information in the bed file and the second one has cds 
#*        information only.
#*					
#*
getAnnotationFromBedFile <- function(bed_file) {
    gene_ranges <- import(bed_file, format="BED");
    if(is.null(gene_ranges$thick))
        stop("cds ranges are not found in bed file.");
		
    cds_ranges <- GRanges(
        seqnames=seqnames(gene_ranges),
        ranges=IRanges::IRanges(start=start(gene_ranges$thick),
            end=end(gene_ranges$thick)),
        strand=strand(gene_ranges));

    names(gene_ranges) <- gene_ranges$name;
    names(cds_ranges)  <- gene_ranges$name;
	
    return (list(gene_ranges,cds_ranges));
}


#* ===========================================================================
#*
#*    3.  getAlignmentsFromBamFile()
#*
#*    Load all alignments from bam file 
#*
#*    Augument: 
#*
#*        bam_file: character vector, a bam file name (and path).
#*
#*    Return:
#*
#*        Granges object containing all reads in bam file. Keep the reads
#*        unfiltered to make this function generic.
#*
getAlignmentsFromBamFile <- function(bam_file) {
    reads <- scanBam(bam_file);

    alignments <- GRanges(
        seqnames=reads[[1]]$rname,
        ranges=IRanges::IRanges(start=reads[[1]]$pos, width=reads[[1]]$qwidth),
        strand=reads[[1]]$strand
    )
	
    rm(reads);
	
    return (alignments);
}


#*  ===========================================================================
#*
#*    4.  countingMetagenePosition()
#*
#*    Find the positions relative to start and stop positions for all 
#*    alignments. Note: genomic positions of alignments and gene/cds 
#*    are all 1-based on forward strand only. 
#*
#*    Arguments:
#*	
#*        alignments:   GRanges object of short reads.
#*        annotation:   GRanges object of cds annotations.
#*        shift_start:  Positive integer, shift this distance. 
#*                      to left from start position of cds.
#*        shift_end:    Positive integer, shift this distance. 
#*                      to right from end position of cds.
#*
#*    Return:
#*
#*        A data frame with 3 columns for distance to cds start, distance
#*        to cds end, and read length for each reads in aligments.	
#*
countingMetagenePosition <- function(alignments, annotation, 
    shift_start=0, shift_end=0)	{
    gene_ranges <- annotation[[1]];
    cds_ranges  <- annotation[[2]];
    meta_cds <- gene_ranges;
	
    plus <- which(strand(meta_cds) == "+");
    start(meta_cds[plus]) <- pmax(0, start(meta_cds[plus]) - shift_start);
    end(meta_cds[plus])   <- end(meta_cds[plus]) + shift_end;
	
    minus <- which(strand(meta_cds) == "-");
    start(meta_cds[minus]) <- pmax(0, start(meta_cds[minus]) - shift_end);
    end(meta_cds[minus])   <- end(meta_cds[minus]) + shift_start;
	
    #*  Reads on forward (plus) strand
    #*  ------------------------------------------------------------
    plus_aligns   <- alignments[strand(alignments) == "+"];
    plus_overlap  <- as.matrix(findOverlaps(plus_aligns, meta_cds));
    plus_to_start <- start(plus_aligns[plus_overlap[,1]]) - 
        start(cds_ranges[plus_overlap[,2]]);
    plus_to_end <- start(plus_aligns[plus_overlap[,1]]) - 
        end(cds_ranges[plus_overlap[,2]]);
    plus_width  <- width(plus_aligns[plus_overlap[,1]]);

    #*  Reads on reverse (minus) strand
    #*  ------------------------------------------------------------
    minus_aligns   <- alignments[strand(alignments) == "-"];
    minus_overlap  <- as.matrix(findOverlaps(minus_aligns, meta_cds));
    minus_to_start <- end(cds_ranges[minus_overlap[,2]]) - 
        end(minus_aligns[minus_overlap[,1]]);
    minus_to_end   <-  start(cds_ranges[minus_overlap[,2]]) - 
        end(minus_aligns[minus_overlap[,1]]);
    minus_width <- width(minus_aligns[minus_overlap[,1]]);
	
    positionProfile <- data.frame(
        to_start=c(plus_to_start, minus_to_start), 
        to_end=c(plus_to_end, minus_to_end), 
        read_len=c(plus_width,minus_width) );

    return (positionProfile);
}


#*  ===========================================================================
#*
#*    5.  countingFrames()
#*
#*    Find frames for all alignments.
#*		
#*    Arguments:
#*	
#*        alignments:  GRanges object of short reads.
#*        cds_ranges:  GRanges object of cds annotations.
#*        inset_5:     positive integer, shift distance after start 
#*                     position of cds to avoid start codon.
#*        inset_3:     positive integer, shift distance before end
#*                     position of cds to avoid stop codon.
#*
#*    Return:
#*
#*        A data frame of two columns for frame and read length.
#*
#*

countingFrames <- function(alignments, cds_ranges, inset_5=34, inset_3=31) {

    #*	Apply insets to avoid start and stop codons
    #*	-------------------------------------------------
    meta_ranges <- cds_ranges; 
	
    plus_cds <- which(strand(meta_ranges) == "+")
    start(meta_ranges[plus_cds]) <- start(meta_ranges[plus_cds]) + inset_5; 
    end(meta_ranges[plus_cds])   <- pmax(start(meta_ranges[plus_cds]), 
        end(meta_ranges[plus_cds]) - inset_3);
	
    minus_cds <- which(strand(meta_ranges) == "-")
    start(meta_ranges[minus_cds]) <- start(meta_ranges[minus_cds]) + inset_3; 
    end(meta_ranges[minus_cds])   <- pmax(start(meta_ranges[minus_cds]), 
        end(meta_ranges[minus_cds]) - inset_5);

    #*  Forward strand reads
    #*  -------------------------------------------------
    plus_aligns <- alignments[strand(alignments) == "+"];
    plus_width <- width(plus_aligns);	
    end(plus_aligns) <- start(plus_aligns);
	
    plus_overlaps <- findOverlaps(plus_aligns, meta_ranges);
    plus_overlaps <- as.matrix(plus_overlaps);
	
    distance <- start(plus_aligns[plus_overlaps[,1]]) - 
            start(cds_ranges[plus_overlaps[,2]]);
    plus_frame <- factor(distance %% 3, levels= 0:2);
    plus_width <- plus_width[plus_overlaps[,1]];	

    #*	reverse strand reads
    #*	---------------------------------------------------	
    minus_aligns <- alignments[strand(alignments) == "-"];
    minus_width <- width(minus_aligns);		
    start(minus_aligns) <- end(minus_aligns);
	
    minus_overlaps <- findOverlaps(minus_aligns, meta_ranges);
    minus_overlaps <- as.matrix(minus_overlaps);
    distance <- end(cds_ranges[minus_overlaps[,2]]) - 
        end(minus_aligns[minus_overlaps[,1]]); 		
    minus_frame <- factor(distance %% 3, levels = 0:2);
    minus_width <- minus_width[minus_overlaps[,1]];	
	
    frame_profile <- data.frame(frame=c(plus_frame, minus_frame), 
        read_len=c(plus_width, minus_width));
	
    return (frame_profile);
}


#*  ===========================================================================
#*
#*    6.  getPositionProfileTable()
#*
#*    Convert position profile table to metagene table.
#*
#*    Arguments:
#*
#*        posProfile:  A Data frame with to_start, to_end, and read_len.
#*        at_which:    A positive integer of 1 or 2.
#*        meta_start:  Integer, start position of metagene
#*        meta_end:    Integer, stop position of metagene
#*
#*    Return:
#*
#*        A numeric matrix with rows for position and columns for 
#*        read length.
#*

getPositionProfileTable <- function(positionProfile=NULL, 
    at_which=1, meta_start=-100, meta_end=100) {

    if(is.null(positionProfile)) stop("Missing profile table.");
    if(!at_which %in% 1:2) stop("at_which must be 1 or 2");
	
    read_len  <- unique(positionProfile[,3]);
    read_len  <- sort(read_len);
    positions <- meta_start:meta_end;
	
    posTable <- initializeMetageneTable(meta_start, meta_end,
        min(read_len), max(read_len));
	
    for(a_row in 1:length(positions)){
        rows <- which(positionProfile[, at_which] == positions[a_row]);
        if(length(rows) > 0) {
            lengths <- positionProfile[rows, 3];
            for(a_col in 1:length(read_len)){
                total <- length(which(lengths == read_len[a_col]));
                posTable[a_row, a_col] <- total;
            }
        }
    }
	
    at_pos_table <- data.frame(total=rowSums(posTable, na.rm=TRUE), posTable);
    colnames(at_pos_table) <- c("total", colnames(posTable));
	
    return(at_pos_table);
}


#*  ===========================================================================
#*
#*    7.  initializeMetageneTable()
#*
#*    Initialize an empty matrix with rows for base positions of metagene 
#*    and columns for read length.
#*
#*    Arguments:
#*
#*        meta_start:  Negative integer, distance before start codon.
#*        meta_end:    Positive integer, distance after end codon.
#*        min_len:     Positive integer, minimum length of alignments.
#*        max_len:     Positive integer, maximum length of alignments.
#*
#*    Return:  A numeric matrix with all zeros.
#*
#*    Note: the total count column is not included and will be added later.
#*
initializeMetageneTable <- function(meta_start=-100, meta_end=100, 
    min_len=25, max_len=34) {

    #*  Augument validation
    #*  --------------------------------------------------
    if(meta_start > 0) stop("metagene start must be negative.");
    if(meta_end   < 0) stop("meta gene end must be positive.");
    if(min_len < 0) stop("min_len must be positive.");
    if(max_len < 0) stop("max_len must be positive.");
    if(max_len < min_len) stop("min_len smaller than max_len.");

    #*  Total row and column numbers
    #*  --------------------------------------------------
    total_rows <- length(meta_start:meta_end);
    total_cols <- length(min_len:max_len);
    profileTable <- matrix(rep(0, total_rows*total_cols), 
        ncol=total_cols);
	
    #*  Coulmn headers
    #*  --------------------------------------------------
    colnames(profileTable) <- min_len:max_len;
    rownames(profileTable) <- meta_start:meta_end;

    return (profileTable);
}


#*  ===========================================================================
#*
#*    8.  getFrameTable()
#*	
#*    Convert frameProfile table to frame table with rows for each read
#*    length and columns for total number and fraction of each frame.
#*
#*    Argument:
#*		
#*        frameProfile:	A data frame with two columns for frames and read_len.
#*
#*    Return:
#*
#*        A data frame with 7 columns for read_len, fraction, frame0  counts,  
#*        frame1 counts, frame2 counts and fraction of the three frames.
#*
getFrameTable <- function(frameProfile) {

    read_len <- unique(frameProfile[,2]);
    read_len <- sort(read_len);
	
    zeros <- frameProfile[frameProfile[,1] == 1, 2];
    ones  <- frameProfile[frameProfile[,1] == 2, 2];
    twos  <- frameProfile[frameProfile[,1] == 3, 2];
	
    frame0 <- rep(0, length(read_len));	
    frame1 <- rep(0, length(read_len));
    frame2 <- rep(0, length(read_len));
	
    for(a_len in 1:length(read_len)) {
        frame0[a_len] <- length(which(zeros == read_len[a_len]));
        frame1[a_len] <- length(which(ones == read_len[a_len]));
        frame2[a_len] <- length(which(twos == read_len[a_len]));
	}

    total <- nrow(frameProfile);	
    len_total <- rowSums(cbind(frame0, frame1, frame2));
    fraction_total  <- signif(len_total/total, digits=4);
    fraction_frame0 <- signif(frame0/len_total, digits=4);
    fraction_frame1 <- signif(frame1/len_total, digits=4);	
    fraction_frame2 <- signif(frame2/len_total, digits=4);

    frame_table <- data.frame(
        fraction=fraction_total, 
        Frame0=frame0,Frame1=frame1,Frame2=frame2,
        Fract0=fraction_frame0,
        Fract1=fraction_frame1,
        Fract2=fraction_frame2, row.names=read_len)
	
    return (frame_table);
}


#*  ===========================================================================
#*
#*    9.  getASiteTable()
#*
#*    Generate a-site table from the three position profile tables
#*
#*    Arguments:
#*
#*        frame_table:    Data frame retured from getFrameTable()
#*        at_stat_table:  Data frame returned from getPositionProfileTable()
#*        at_end_table:   Data frame returned from getPositionProfileTable() 
#*        parameters:     List of length 13, all parameters for framing
#*
#*    Return:
#*
#*        A data frame with columns for:
#*
#*        read length
#*        fraction of each length,
#*        at_start (startPeak)
#*        at_end (endPeak)
#*        info
#*        fraction of of each frame at each length
#*        Asite (bestASite)
#* 
#*
getASiteTable <- function(frame_table, at_start, at_end, parameters) {

    minLenFract <- parameters$minLenFract;
    minStartRange	<- parameters$minStartRange; 
    maxStartRange	<- parameters$maxStartRange;
    minEndRange <- parameters$minEndRange; 	
    maxEndRange <- parameters$maxEndRange;	
    startShift	<- parameters$startShift;	
    endShift	<- parameters$endShift;		

    #*  atStart and atEnd columns
    #*  ----------------------------------------------
    start_peak <- startShift + getReadPeaks(at_start, minStartRange, 
        maxStartRange);
    end_peak <- endShift + getReadPeaks(at_end, minEndRange, 
        maxEndRange);
						
    #*  Info column
    #*  ------------------------------------------------------------
    Info <- log2(3) - rowSums(cbind(frame_table[,5]*log2(frame_table[,5]), 
        frame_table[,6]*log2(frame_table[,6]), 
        frame_table[,7]*log2(frame_table[,7]) ) ) * -1

    #*  ASite column
    #*  ----------------------------------------------------
    start_frame <- (start_peak*-1) %% 3;
    end_frame   <- (end_peak*-1) %% 3;
    best_frame  <- getBestFrame(frame_table);
	
    best_asite <- rep(0, nrow(frame_table));
    for(a_row in 1:length(best_asite)) {
        if(start_peak[a_row] == end_peak[a_row] && 
            best_frame[a_row] == start_frame[a_row]) {
            best_asite[a_row] <- start_peak[a_row];
        } else if(end_frame[a_row] == best_frame[a_row]) {
            best_asite[a_row] <- end_peak[a_row];
        } else if (start_frame[a_row] == best_frame[a_row]) {
            best_asite[a_row] <- start_peak[a_row];
        } else if (start_peak[a_row] == end_peak[a_row]) {
            best_asite[a_row] <- paste0(start_peak[a_row],"???");
        } else {
            best_asite[a_row] <- "***";
        }
    }
	
    a_site_table <- data.frame(
        Fract=frame_table[,1],
        AtStart=start_peak, AtEnd=end_peak,
        Info=Info, Frame0=frame_table[,5],
        Frame1=frame_table[,6],
        Frame2=frame_table[,7],
        ASite=best_asite, 
        row.names=rownames(frame_table));

    a_site_table <- a_site_table[a_site_table[,1] > minLenFract, ];
	
    return (a_site_table);
}


#*  ===========================================================================
#*
#*    10.	getReadPeaks()
#*
#*    Find count peaks for each read length from metagene length table.
#*			
#*    Arguments:
#*
#*        metagene:   A data frame returned from getPositionProfileTable()  
#*                    with columns for each length of reads and rows for 
#*                    counts in each position of metagene.
#*        min_range:  Positive integer, minimum index of the row to find peak.
#*        max_range:  Positive integer, maximum index of the row to find peak.
#*
#*    Return:
#*
#*        A integet vector, distance from read start for each read length.
#*
#*
getReadPeaks <- function(metagene, min_range, max_range) {

    #*  working on subset (rows min_range:max_range)
    #*  ============================================
    row_ranges <- min_range:max_range;
    rows <- which(as.numeric(rownames(metagene)) %in% row_ranges);
    count_range <- as.matrix(metagene[rows, 2:ncol(metagene)]);
	
    peaks <- rep(0, ncol(count_range));
    max_count <- apply(count_range, 2, max);
	
    for(a_col in 1:ncol(count_range)){
        the_index <- which(count_range[,a_col] == max_count[a_col]);
        if(length(the_index) > 1) 
            the_index <- the_index[length(the_index)];
        peaks[a_col] <- the_index;		
    }
    distance <- (row_ranges*-1)[peaks];
	
    return (distance);
}


#*  ===========================================================================
#*
#*    11.  getBestFrame()
#*
#*    Find the best frame from Frame 0~2 for each read length
#*
#*    Argument:
#*
#*        An data frame with counts of 3 frames for each read length
#*
#*    Return:
#*
#*        Positive integer vector, frame number for each read length,
#*
#*
getBestFrame <- function(frame_table) {

    frame_counts <- frame_table[,2:4]
    best_frame <- rep(0, nrow(frame_counts));
    for(a_len in 1:length(best_frame)) {
        the_index <- which(frame_counts[a_len,] == max(frame_counts[a_len,]));
        best_frame[a_len] <- the_index - 1;
    }

    return (best_frame);
}


#* ===========================================================================
#*
#*    12.  validateParameters()
#*
#*    Validate parameters in case user edited the default one
#*
#*    Argument: parameters, list of numeric variables
#*
#*    Return: None
#*
validateParameters <- function(parameters) {

    if(parameters$shift_start < 0) stop("shift_start must be >= 0.");
    if(parameters$shift_end < 0 )  stop("shift_end must be >= 0.");
    if(parameters$inset_5 < 0)     stop("inset_5 must be >= 0.");
    if(parameters$inset_3 < 0)     stop("inset_3 must be >= 0.");
    if(parameters$min_len < 20)    stop("min_len must be >= 20");
    if(parameters$max_len < 30)    stop("max_len must be >= 30");
    if(parameters$startShift < 0)  stop("startShift must be > 0.");
    if(parameters$endShift > 0)    stop("endShift must be < 0");
	
    if(parameters$minStartRange > 0)
        stop("minStartRange must be less than 0.");
    if(parameters$maxStartRange > 0) 
        stop("maxStartRange must be less than 0.");
    if(parameters$minEndRange > 0) 
        stop("minEndRange must be less than 0.");
    if(parameters$maxEndRange > 0) 
        stop("maxEndRange must be less than 0.");

    if(parameters$minLenFract > 1 || parameters$minLenFract < 0) 
        stop("minLenFract must be between 0 and 1.");	
}


#*	===========================================================================
#*
#*    13.  getMetageneFrames()
#*
#*    Summrize reading frames for each qualified reads based on metagene 
#*    positions (positions of read start relative to cds start)
#*
#*    Arguments:
#*
#*        posProfile:     A data frame with three columns for to_start,  
#*                        to_end, and read_len
#*
#*        bam_file_name:  Character vector, name of bam file from which the  
#*                        reads are scaned for metagene positions. Used for
#*                        output file generation.
#*
getMetageneFrames <- function(posProfile, bam_file_name) {

    message("Counting frames by distance to cds start.");
    frames_by_start <- summerizeMetageneFrames(posProfile, 1);
    output_file <- sub("bam", "meta_start_frames.txt", basename(bam_file_name));
    write.table(frames_by_start, file=output_file, sep="\t", quote=FALSE,
        row.names=FALSE, col.names=TRUE);

    message("Counting frames by distance to cds stop.");			
    frames_by_stop <- summerizeMetageneFrames(posProfile, 2);
    output_file <- sub("bam", "meta_stop_frames.txt", basename(bam_file_name)); 
    write.table(frames_by_stop, file=output_file, sep="\t", quote=FALSE,
        row.names=FALSE, col.names=TRUE);
}


#*	===========================================================================
#*
#*    14.  summerizeMetageneFrames()
#*
#*    Summrize reading frames for each qualified reads based on metagene 
#*    positions (positions of read start relative to cds start or end)
#*
#*    Arguments:
#*
#*        posProfile:  A data frame with four columns for to_start,  
#*                     to_end, read_len, and frames for each read.
#*
#*        pos_col:     Integer, column numberfor target metagene positions.
#*
#*
#*    Return:
#*
#*        A matrix of 3 columns for each reading frame and rows for each
#*        metagene position
#*	
summerizeMetageneFrames <- function(posProfile, pos_col=1) {

    if(pos_col < 1 || pos_col > 2)
        stop("Incorrect column for metagene positions.");

    positions <- posProfile[, pos_col];
    pos_frames <- positions %% 3;
	
    meta_frames <- initializeMetageneFrameTable(positions);
    meta_position <- as.character(meta_frames$metaPosition);
	
    for(a_frame in 0:2) {
        frame_rows <- which(pos_frames == a_frame);
        counts <- as.data.frame(table(positions[frame_rows]), 
            stringsAsFactors=FALSE);

        rows <- match(counts$Var1, meta_position);
        meta_frames$counts[rows] <- counts$Freq;
    }
	
    if(pos_col == 1) {
        colnames(meta_frames)[3] <- "to_start";
    } else { colnames(meta_frames)[3] <- "to_stop"; }
	
    return (meta_frames);
}




#*  ===========================================================================
#*
#*    15.  initializeMetageneFrameTable()
#*
#*    Initialized a empty matirx(all zeros) with 3 columns and rows of
#*    number of all unique positions (argumnents)
#*
#*    Argument:
#*
#*        position_set:  Integer vector, unique position set of all metagene
#*                       positions (position of reads relative to cds start 
#*                       or cds stop positions)
#*
#*    Return:
#*
#*        A data frame of 3 columns for number of frames, frames, and
#*        metagene positions.
#*
initializeMetageneFrameTable <- function(position_set) {

    pos_set <- unique(position_set);	
    pos_set <- pos_set[order(pos_set)];
	
    #*  each unique position can belong to one frame only 
    #*  =================================================
    frame_table <- data.frame(
        counts=rep(0, length(pos_set)),
        frames=factor(pos_set %% 3, levels=0:2),		
        metaPosition=pos_set
    );
    return (frame_table);
}


#*
#*  End of RiboProfiling_FpFraming.R
#*  ________________________________________________________________________
#*  <fp Framing><fp Framing><fp Framing><fp Framing><fp Framing><fp Framing>


















