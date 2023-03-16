#*
#*    File:  RiboProfiling_Wiggle_Track.R
#*
#*    Functions to generate wiggle files from bam file.
#*
#*    Base on the source code from Dr. Nocholas Ingoliad for ribosomal  
#*    footprinting analsis (https://github.com/ingolia-lab/RiboSeq): 
#*
#*    The key difference of this tool with other software are counting  
#*    coverage for forward and reverse strand separately and the coverage 
#*    is on the A-site only.
#*
#*    Original command man page:
#*
#*    wiggle-track [OPTIONS] <BAM>
#*    -o OUTFILE    --output=OUTFILE     Output filename
#*    -a ASITEFILE  --asite=ASITEFILE    A site offsets filename
#*    -c            --coverage           Total read coverage
#*    -q QNORM      --qnorm=QNORM        Multiplicative scaling factor
#*    -x SIZEFILE   --chrsizes=SIZEFILE  Chromosome size output file
#*
#*
#*    Function implemented:
#*
##    1.  getWiggleCounts()
##    2.  getChromosomeSizesFromBam()
##    3.  writeChromosomeSizesToFile()
##    4.  readChromosomeSizesFromFile()
##    5.  countAtAsite()
##    6.  writeWiggleFiles()
##    7.  normalizeWiggleCounts()
##    8.  sortTableByChromosomeNames()
##    9.  readAsiteFromFile()
##    10. getBigWiggleFileForTranscripts()
##    11. filterOutAlignments()
##    12. getTranscriptReadsOnly()
##    13. mapASite()
##    14. getBigWigFilesByAsite()
##
##    Last revised on October 04, 2022
#*  ______________________________________________________________________
#*  <Wiggle_Track><Wiggle_Track><Wiggle_Track><Wiggle_Track><Wiggle_Track>
#* 




###############################################################################
##
##    1.  getWiggleCounts()
##
##    The main function to get counts on A-site for wiggle tracks. 
##
##    Arguments:
##
##        bam_file:       chraracter vector, name of bam file (with path) 
##        asite_table:    data frame of column 1 for a site and rownames
##                        for read length. A file name is also accepted.
##        annot_bed_file: chraracter vector, name of bed file (with path)
##
##    Return:
##
##        List of list with each sub-list is wiggle counts for 
##        forward and reverse strand of one chromomes.
##

getWiggleCounts <- function(bam_file, asite_table, annot_bed_file=NULL) {

    #*  Get chromosome names and lengths
    #*  ------------------------------------------------
    chrom_length <- getChromosomeSizesFromBam(bam_file);

    #*  Read in aiste from file if have not done so
    #*  ------------------------------------------------
    if(is.character(asite_table)) 
        asite_table <- readAsiteFromFile(asite_table);

    #*	Read in all reads from bam file and filter out
    #*	un-wanted ones)
    #*	------------------------------------------------
    message("Read aligments from ", bam_file);
    alignments <- getAlignmentsFromBamFile(bam_file);
    alignments <- filterOutAlignments(alignments, 
        asite_table, annot_bed_file);

    #*	Count coverage at a-site for each chromosome
    #*  ------------------------------------------------
    all_counts <- list();
    for(a_chrom in 1:nrow(chrom_length))
    {
        chrom_name <- as.character(chrom_length[a_chrom, 1]);
        chrom_len <- as.numeric(chrom_length[a_chrom, 2]);		
        message(paste("count coverage for", chrom_name));

        chrom_index <- which(seqnames(alignments) == chrom_name);
        if(length(chrom_index) == 0) next;
		
        all_counts[[chrom_name]] <- countAtAsite(alignments[chrom_index], 
            chrom_len, asite_table);	
    }

    return (all_counts);
}


###############################################################################
##
##    2.  getChromosomeSizesFromBam()
##
##    Extract chromosome length from bam file header. It is supposed
##    here all bam files have header lines with chromosome length. The
##    chromosome names and order are as is to leave it for downstream
##    processing.
##
##    Argument:
##
##        bam_file, character vector, bam file name (and path).
##
##    Return:
##
##        A data frame of columns 2 for chromosome names and lengths.
##

getChromosomeSizesFromBam <- function(bam_file)
{
    bam_headers <- scanBamHeader(bam_file, what=c("targets"));
    chromSizes  <- bam_headers[[1]]$targets;

    chromSizes <- data.frame(
        chromosome=names(chromSizes),
        chromLength=as.numeric(chromSizes)
    );

    return (chromSizes);
}

writeChromosomeSizesToFile <- function(chromSizes, file_name)
{
    write.table(chromSizes, file=file_name, sep="\t",
        quote=FALSE, row.names=FALSE, col.names=TRUE);
}
readChromosomeSizesFromFile <- function(file_name)
{
    chromSizes <- read.table(file_name, header=TRUE,
        sep="\t", quote="");

    return (chromSizes);
}


###############################################################################
##
##    5.  countAtAsite()
##
##    Count read number overlap at A-site for transcripts on one chromosome. 
##    The overlap is stranded.
##
##    Arguments:
##
##        alignments:   GRange object for short reads on one chromosome.
##        chrom_len:    Positive integer, length of the chromosome.
##        asite_table:  data frame of one column for A-sites and rownames
##                      for qualified read length.
##
##    Return:	
##
##        A list of 2 integer vectors for counts on forward and reverse
##        strand at each base pair position
##
countAtAsite <- function(alignments, chrom_len, asite_table)
{
    #*	Coverage at base level on forward strand
    #*	-------------------------------------------------------
    fwd_aligns <- alignments[which(strand(alignments) == "+")];
    fwd_aligns <- mapASite(fwd_aligns, asite_table);

    fwd_counts <- rep(0, chrom_len);
    locations <- start(fwd_aligns);
    counts <- table(locations);
    count_index <- as.numeric(names(counts));
    fwd_counts[count_index] <- counts;

    #*	Coverage at base level on reverse strand
    #*	-------------------------------------------------------	
    rev_aligns <- alignments[which(strand(alignments) == "-")];
    rev_aligns <- mapASite(rev_aligns, asite_table); 
	
    rev_counts <- rep(0, chrom_len);
    locations <- start(rev_aligns);
    counts <- table(locations);
    count_index <- as.numeric(names(counts));
    rev_counts[count_index] <- counts;
	
    count_list <- list(forward=fwd_counts, reverse=rev_counts);
    return (count_list);
}



###############################################################################
##
##    6.  writeWiggleFiles()
##
##    Write wiggle counts to wiggle track file
##
##    Arguments:
##
##        all_counts:    list of list, each list element has two numeric
##                       vectors for read counts forward and reverse 
##                       strand at base pair level.
##        bam_file:      character vector, name of bam file from which 
##                       the wiggle track data is generated.
##        is.normalized: Logic, if the counts are normalized.
##
##    Return:	none. Write files only
##
##
writeWiggleFiles <- function(all_counts, bam_file, is.normalized=FALSE)
{
    if(is.normalized) {
        forward_file <- sub("bam$", "fwd_normed.wig", bam_file);
        reverse_file <- sub("bam$", "rev_normed.wig", bam_file);
    } else {
        forward_file <- sub("bam$", "forward.wig", bam_file);
        reverse_file <- sub("bam$", "reverse.wig", bam_file);
    }

    if(file.exists(forward_file)) file.remove(forward_file);
    if(file.exists(reverse_file)) file.remove(reverse_file);

    file.create(forward_file);
    file.create(reverse_file);

    for(a_chrom in 1:length(all_counts)) {
        chrom_name <- names(all_counts)[a_chrom];
        header <- paste0("variableStep chrom=", chrom_name, " span=1");

        fwd_count <- all_counts[[chrom_name]]$forward;
        fwd_count <- data.frame(1:length(fwd_count), fwd_count);
        fwd_count <- fwd_count[which(fwd_count[,2] > 0),];
        fwd_count <- paste(fwd_count[,1], fwd_count[,2]);

        write(header, file=forward_file, append=TRUE);
        write(fwd_count, file=forward_file, append=TRUE);

        rev_count <- all_counts[[chrom_name]]$reverse;
        rev_count <- data.frame(1:length(rev_count), rev_count);
        rev_count <- rev_count[which(rev_count[,2] > 0),];
        rev_count <- paste(rev_count[,1], rev_count[,2]);

        write(header, file=reverse_file, append=TRUE);
        write(rev_count, file=reverse_file, append=TRUE);
    }
}


###############################################################################
##
##    7.  normalizeWiggleCounts()
##
##    Normalize the wiggle track data to a defined total counts
##
##    Arguments:
##
##        all_counts:       list of list, each list element has two 
##                          numeric vectors for read counts forward  
##                          and reverse strand at base pair level.
##        normalize_factor: positive numeric, scaling factor.
##
##    Return:
##
##        list of list, all normalized counts 
##
##
normalizeWiggleCounts <- function(all_counts, normalize_factor=NULL) {

    #*	If normlaize factor is not provided
    #*	---------------------------------------
    if(is.null(normalize_factor)) 
    {
        ## size factor based on 1 million reads with 
        ## read length of 100
        scale_to <- 100000000;

        wig_total <- 0;
        for(a_chrom in 1:length(all_counts))
        {
            wig_total <- wig_total + sum(all_counts[[a_chrom]]$forward) + 
                sum(all_counts[[a_chrom]]$reverse);
        }
        norm_f <- scale_to/wig_total;
    } else { norm_f <- normalize_factor; }

    for(chrom in 1:length(all_counts))
    {
        normed_fwd <- all_counts[[chrom]]$forward*norm_f;
        all_counts[[chrom]]$forward <- round(normed_fwd, digits=2);
		
        normed_rev <- all_counts[[chrom]]$reverse*norm_f;
        all_counts[[chrom]]$reverse <- round(normed_rev, digits=2);
    }

    return (all_counts);
}


###############################################################################
##
##    8.  sortTableByChromosomeNames
##
##    Sort a table by chromosome names (either Arabic or Roman numbers)
##
##    Arguments:
##
##        chrom_info: data frame or matrix with column of chromosome
##						names.
##        name_col:   positive integer, column number of chromosome names
##        type:       character vector, type of chromosome numbers, 
##                    either "digit" or "roman".
##
##    Return:
##
##        Data frame or matrix, sorted input table
##
##
sortTableByChromosomeNames <- function(chrom_info, name_col=1, type="digit") {

    digit_char <- as.character(1:50);
    roman_char <- as.character(as.roman(1:50));

    chrom_name <- as.character(chrom_info[,name_col]);
    if (length(chrom_name) != length(unique(chrom_name))) 
        stop("Chromosome names is not an unique set.\n")

    #*	Remove "chr" from chromosome names, if it exists
    #*	------------------------------------------------------
    with_prefix <- length(grep("^chr", chrom_name));
    if (with_prefix > 0) {
        if (with_prefix != length(chrom_name)) {
            stop("Not all chromosome name have prefix 'chr'");
        }
        else { chrom_name <- gsub("^chr", "", chrom_name);}
    }

    #*	Only sort the numbered chromosomes
    #*	------------------------------------------------------
    type <- tolower(type);
    if(type == "digit") {
        digit_rows <- which(chrom_name %in% digit_char);
        chrom_digit <- as.numeric(chrom_name[digit_rows]);
        chrom_table <- chrom_info[digit_rows, ];
        chrom_table <- chrom_table[order(chrom_digit),];
		
        chrom_other <- chrom_info[-digit_rows, ];
        chrom_table <- rbind(chrom_table, chrom_other);

    } else if(type == "roman") {
        roman_rows <- which(chrom_name %in% roman_char);
        chrom_digit <- as.integer(as.roman(chrom_name[roman_rows]));		
        chrom_table <- chrom_info[roman_rows, ];
        chrom_table <- chrom_table[order(chrom_digit),];

        chrom_other <- chrom_info[-roman_rows, ];
        chrom_table <- rbind(chrom_table, chrom_other);
    } else { stop("Incorrect chromosome name type."); }

    return (chrom_table);
}


###############################################################################
##
##    9.  readASiteFromFile()
##
##    Read a-site table generated by FpFraming(). Since the original code 
##    from Ingolia's source file uses colname names ('asite') to refer
##    the a-site values for each read length, the exact same column name
##    must be set. Otherwise, getAsiteProfile() will return all 0's.
## 
##    Argument:
## 
##        Character vector, name (and path) of a-site file generated in the
##        step of FpFraming().
##
##    Return:
##
##        A data frame of 1 columns to hold a-site for each read length which
##        are represented by row names.
##

readAsiteFromFile <- function(a_site_file) {
    asite_table <- read.table(a_site_file, header=TRUE, 
        sep="\t", quote="");
    asites <- data.frame(asite=asite_table[,2], 
        row.names=rownames(asite_table));

    return (asites);
}


###############################################################################
##
##    10. getBigWiggleFileForTranscripts()
##
##    Simply generated bigWig file with rtracklayer and GenomicAlignment
##    pacakges.
##
##    Arguments:
##
##        bam_file:    character vector, name of bam file from which 
##                     the wiggle track data is generated.
##        transcripts: GRange object for a transcript.
##
##    Return: None
##
##    Last revised on January 7, 2019
##
getBigWiggleFileForTranscripts <- function(bam_file, transcripts) {
    alignments <- getAllAlignments(bam_file, transcripts)

    forward_gr <- alignments[which(strand(alignments) == "+")]
    forward_cov <- coverage(forward_gr);
    forward_file <- sub("bam$", "_forward.bigwig", basename(bam_file));
    export.bw(forward_cov, forward_file);

    reverse_gr <- alignments[which(strand(alignments) == "-")];
    reverse_cov <- coverage(reverse_gr);
    reverse_file <- sub("bam$", "_reverse.bigwig", basename(bam_file));
    export.bw(reverse_cov, reverse_file);
}


###############################################################################
##
##    11.  filterOutAlignments()
##
##    Filter out alignments by read length and annotation regions.
##
##    Arguments:
##
##        alignments:     GRange object of all alignments
##        asite_table:    A data frame of 1 columns to hold a-site for  
##                        each read length which are represented by 
##                        row names.
##        annot_bed_file: character vector, name (and path) of annotation    
##                        file in bed format.
##
##    Return: 
##
##        GRange object of filtered alignments.
##
##
##    Last revised on January 7, 2019
##

filterOutAlignments <- function(alignments, asite_table, annot_bed_file) {
    #*	Filter out the reads by length
    #*	------------------------------------------------

    read_len <- as.numeric(rownames(asite_table));	
    align_width <- width(alignments);

    bad_reads <- c(which(align_width < min(read_len)), 
                   which(align_width > max(read_len)));
    if(length(bad_reads) > 0)
        alignments <- alignments[-bad_reads];

    #*	If only transcript regions are interested
    #*	---------------------------------------------------------
    if(!is.null(annot_bed_file))
        alignments <- getTranscriptReadsOnly(alignments, annot_bed_file);

    return (alignments);
}


###############################################################################
##
##    12.  getTranscriptReadsOnly()
##
##    Filter out alignments by annotation regions.
##
##    Arguments:
##
##        alignments:     GRange object of all alignments
##        annot_bed_file: Character vector, name (and path) of annotation    
##                        file in bed format.
##
##    Return:  
##
##        GRange object of all alignments.
##
##    Last revised on January 7, 2019
##
getTranscriptReadsOnly <- function(alignments, annot_bed_file) {

    annot_bed <- import(annot_bed_file, "BED");
    names(annot_bed) <- annot_bed$name;

    #*  Check out if the annotation file is correct
    #*  -------------------------------------------
    bam_chrom <- unique(seqnames(alignments));
    bed_chrom <- unique(seqnames(annot_bed));

    if(length(which(bam_chrom %in% bed_chrom)) != length(bam_chrom))
        message("Some chromosome is missing in annotation file.")

    overlaps <- findOverlaps(annot_bed, alignments);
    overlaps <- as.matrix(overlaps);
    if(nrow(overlaps) == 0) stop("No overlap found.");

    transcript_reads <- alignments[as.numeric(overlaps[,2])];

    return (transcript_reads);
}


###############################################################################
##
##    13. mapASite()
##
##    Convert alignments to one base width GRanges based on a-site.
##
##    Arguments:
##
##        alignments:  GRange object of all alignments
##        asite_table: A data frame of 1 columns to hold a-site for each  
##        read length which are represented by row names.
##
##    Return:
##
##        GRange object of alignments overlapped with transcript
##        regions only.
##
##    Last revised on January 7, 2019
##

mapASite <- function(alignments, asite_table)
{
    align_width <- width(alignments);
    align_strand <- as.character(unique(strand(alignments)));

    read_len <- as.numeric(rownames(asite_table));
    a_site <- as.numeric(asite_table[,"asite"]);

    for(a_row in seq_along(read_len))
    {	
        A_offset <- a_site[a_row];
        alns_len <- read_len[a_row];	

        if(align_strand == "+") {
            align_start <- A_offset + 1;
        } else {
            align_start <- alns_len-A_offset;
        }

        read_index <- which(width(alignments) == alns_len);
        alignments[read_index] <- narrow(alignments[read_index], 
            start=align_start, width=1);
    }

    return (alignments);
}


###############################################################################
##
##    14. getBigWigFilesByAsite()
##
##    Generate BigWig file from bam file with a-site table.
##
##    Arguments:
##
##        bam_file:     character vector, name (and path) of bam file.
##        bam_type:     character vector, either "Riboseq" or "RNASeq". 
##        asite_file:   character vector, name (and path) of a-site file.
##        norm_method:  character vector, type of normalization.
##        scale_to:     Positive integer, scale total counts to this level.
##                      Default: 10000000 (10 million reads)
##
##    Return: None. Result will be write to file.
##
##    Last modified on December 05, 2022
##

getBigWigFilesForAsite <- function(bam_file=NULL, bam_type="Riboseq", 
    asite_file=NULL, norm_method=c("raw", "scale", "cpm"), 
    scale_to=1000000) {

    if(is.null(bam_file) | is.null(asite_file))
        stop("Missing input file name(s).");

    bam_type <- tolower(bam_type);
    if(bam_type == "riboseq") {
        asite <- readAsiteFromFile(asite_file);
    } else { asite <- getDefaultRNASeqAsites();}

    read_lengths <- as.numeric(rownames(asite))
    min_len <- min(read_lengths);
    max_len <- max(read_lengths);	

    alignments <- getAlignmentsFromBamFile(bam_file);
    align_width <- width(alignments);
    to_keep <- which(align_width >=min_len & align_width <= max_len);
    alignments <- alignments[to_keep];

    fwd_aligns <- alignments[which(strand(alignments) == "+")];
    fwd_aligns <- mapASite(fwd_aligns, asite);

    rev_aligns <- alignments[which(strand(alignments) == "-")];
    rev_aligns <- mapASite(rev_aligns, asite);

    fwd_cov <- coverage(fwd_aligns);
    rev_cov <- coverage(rev_aligns);

    asite_alignments <- c(fwd_aligns, rev_aligns);
    asite_cov <- coverage(asite_alignments);	

    norm_method <- tolower(norm_method);
    if(norm_method == "scale") {
        scale_factor <- scale_to/sum(sum(asite_cov));

        scaled_fwd_cov <- fwd_cov * scale_factor;
        export(scaled_fwd_cov, con=sub("bam", "fwd_strand_scaled.bw", 
            basename(bam_file)), "BW");

        scaled_rev_cov <- rev_cov * scale_factor;
        export(scaled_rev_cov, con=sub("bam", "rev_strand_scaled.bw", 
            basename(bam_file)), "BW");

        scaled_asite_cov <- asite_cov * scale_factor;
        export(scaled_asite_cov, con=sub("bam", "both_strand_scaled.bw", 
            basename(bam_file)), "BW");
    } else if (norm_method == "cpm") {
        norm_factor <- 1000000/length(asite_alignments);

        cpm_fwd_cov <- fwd_cov * norm_factor;
        export(cpm_fwd_cov,  con=sub("bam", "fwd_strand_cpm.bw", 
            basename(bam_file)), "BW");

        cpm_rev_cov <- rev_cov * norm_factor;
        export(cpm_rev_cov,  con=sub("bam", "rev_strand_cpm.bw", 
            basename(bam_file)), "BW");

        cpm_asite_cov <- asite_cov * norm_factor;
        export(cpm_asite_cov,  con=sub("bam", "both_strand_cpm.bw", 
            basename(bam_file)), "BW");
    } else if (norm_method == "raw") {
        export(fwd_cov, con=sub("bam", "fwd_strand_raw_bw", 
            basename(bam_file)), "BW");

        export(rev_cov, con=sub("bam", "rev_strand_raw_bw", 
            basename(bam_file)), "BW");

        export(asite_cov, con=sub("bam", "both_strand_raw_bw", 
            basename(bam_file)), "BW");
    } else { stop("Only 'Raw', 'CPM', 'Scale' method supported."); }

    message("Done for: ", basename(bam_file));
}


##
##  End of RiboProfiling_Wiggle_Track.R
##  ______________________________________________________________________
##  <Wiggle_Track><Wiggle_Track><Wiggle_Track><Wiggle_Track><Wiggle_Track>

