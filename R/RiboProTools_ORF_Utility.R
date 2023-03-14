#*
#*    File:  ORF_Utility.R
#*
#*    Source code for Extracting ORF information from DNA sequence. This  
#*    file could be used as a standlone source file for loading or adding 
#*    to other package source such as RiboPro.R files.
#*
#*    Functions implemented:
#*
#*    1.   getCodonFromSequence()
#*    2.   is.mRNA.sequence()
#*    3.   validateCodons()
#*    4.   getDefaultStartCodon()
#*    5.   getDefaultStopCodon()
#*    6.   getAnnotationFromBigBedFile()
#*    7.   checkBlocksInBed()
#*    8.   getSequenceFromMutipleFastaFiles()
#*    9.   getSequenceFromOneFastaFile()
#*    10.  getUTRInfoDefinedByBed()
#*    11.  checkChromosomeInfo()
#*    12.  getUTRInfoFromFastaFile()
#*    13.  covertToRNASequence()
#*    14.  getReversedSequence()
#*    15.  getComplementarySequence()
#*    16.  is.DNA.sequence()
#*    17.  getCDSInfoDefinedByBed()
#*    18.  writeORFInfoToFile()
#*    19.  getCodonIndex()
#*    20.  setCodonIndexPair()
#*    21.  screenORF()
#*    22.  getORFPositions()
#*    23.  getStartCodonContext()
#*    24.  converToGenePredFormat()
#*
#*    Edited on June 25, 2019 for naming conventions, check possible bugs,
#*    and add function for identifying ORF in CDS sequence.
#*
#*  _____________________________________________________________________
#*  <ORF_Utility>*<ORF_Utility>*<ORF_Utility>*<ORF_Utility>*<ORF_Utility>




#*  ===================================================================
#*		
#*    1.  getCodonFromSequence()
#*
#*    Transfer mRNA sequence to codons. 
#*
#*    Arguments:
#*
#*        mRNA_seq:  Character vector, a fragment of mRNA sequence. 
#*                   Each base must be defined (no "N" allowed).
#*        start_at:  Positive integer, either 1, 2, or 3,  which 
#*                   position where the codons will start
#*
#*    Return:
#*
#*        Character vector, all codons starting from the start_at 
#*        position in the mRNA sequence
#*
#*		Last edited on June 17, 2019
#*
getCodonFromSequence <- function(mRNA_seq, start_at) {

    if(is.mRNA.sequence(mRNA_seq) == FALSE)
        stop("Input must be mRNA sequence.")

    if(is.numeric(start_at) == FALSE)
        stop("The start position must be 1, 2, or 3.")
    if(start_at < 0 || start_at > 3)
        stop("The start position must be 1, 2, or 3.")

    #*  How many codons
    #*  ==========================================

    total_bases <- nchar(mRNA_seq);
    base_one <- seq(start_at, total_bases, 3);

    last_base <- base_one[length(base_one)] + 2;
    if(last_base > total_bases) 
        base_one <- base_one[-length(base_one)];
    
    #*  separate sequences to codons
    #*  ==============================================
   
    codons <- rep("   ", length(base_one));
    for(a_codon in 1: length(base_one))
      codons[a_codon] <- substr(mRNA_seq, 
           base_one[a_codon], base_one[a_codon] + 2);

    return (codons);
}

#*  ===================================================================
#*		
#*    2.  is.mRNA.sequence()
#*
#*    Check out if the sequence is mRNA sequence, meaning the bases
#*    in sequence must be A, U, G, C only.
#*
#*    Argument:
#*    
#*        mRNA_seq: character vector, a fragment of mRNA sequence
#*
#*    Return: 
#*
#*        Logic, TURE if argument is mRNA, otherwise FALSE
#* 
#*    Last edited on June 17, 2019
#*
is.mRNA.sequence <- function(mRNA_seq) {

    if(is.character(mRNA_seq) == FALSE) {
        message("Sequences must be character vector.")
        return (FALSE);
    }
    mRNA_seq <- toupper(mRNA_seq);

    mRNA_bases <- c("A", "U", "G", "C");
    seq_bases <- strsplit(mRNA_seq, split="")[[1]];

    if(sum(seq_bases %in% mRNA_bases) != length(seq_bases)) {
        message("mRNA sequence bases must be A, U, G or C.")
        return (FALSE);
    } else {
        return (TRUE);
    }
}

#*  ===================================================================
#*		
#*    3.  validateCodons()
#*
#*    Check codons. Codons must be held by character vector and each codon
#*    must have 3 bases from "A", "U", "G", "C".
#*
#*    Argument:
#*
#*        codons: character vector which holds one or more codons.
#*
#*    Return: none
#*
#*    Last edited on June 17, 2019
#*
validateCodons <- function(codons) {

    if(is.character(codons) == FALSE)
        stop("Codon must be character vector.")


    codon_base <- c("A", "U", "G", "C");
    for(a_codon in 1:length(codons))
    {
        if(nchar(codons[a_codon]) != 3)
            stop("Each codon must have 3 bases.")

        bases <- strsplit(codons[a_codon], split="")[[1]];
        if(sum(bases %in% codon_base) != 3)
            stop("Codon base must be A, U, G or C.")
    }
}

#*  ===================================================================
#*		
#*    4.  getDefaultStartCodon()
#*
#*    Get default start codon set
#*
#*    Argument:	None.
#*    Return: 	Default codon set.
#*
#*    Last edited on June 17, 2019
#*
getDefaultStartCodon <- function() {

    start_codons <- c("AUG", "AUC", "AUU", "AUA", "AGG",
                     "ACG", "AAG", "CUG", "UUG", "GUG");  
    return (start_codons);
}

#*  ===================================================================
#*		
#*    5.  getDefaultStopCodon()
#*
#*    Get default stop codon sets
#*		
#*    Argument: None.
#*    Return:   Default codon set ("UAG", "UGA", and "UAA").
#*		
#*    Last edited on June 17, 2019
#*		
getDefaultStopCodon <- function() {

    stop_codon <- c("UAG", "UGA", "UAA");
    return (stop_codon);
}

#*  =========================================================================
#*
#*    6.  getAnnotationFromBigBedFile()
#* 
#*    Read bigBed file and check the contents
#*
#*    Arguments: 
#*
#*        bed_file:       character vector, bed file name
#*        has.header:     logic, if the bed file has a header line
#*        sepcial_chrom:  character vector, chromosome names other than 
#*                        digits, roman numbers, and "X", "Y", "M"
#*
#*		Return:	  A data frame with all contents in the bed file
#*
#*    Fields in bigBed file: 
#*    (https://genome.ucsc.edu/FAQ/FAQformat#format1.7)
#*
#*    1).   chrom:      chromosome name of each feature
#*    2).   chromStart: start position of each feature
#*    3).   chromEnd:   end position of each feature
#*    4).   name:       gene name of each feature
#*    5).   score:      used for graphic display only
#*    6).   strand:     chromosome strand of each feature
#*    7).   thickStart: starting position drawn thickly
#*    8).   thickEnd:	 ending position drawn thickly
#*    9).   itemRgb:    used for graphic display only
#*    10).  blockCount: number of blocks of each feature
#*    11).  blockSizes: comma-separated list of block sizes
#*    12).  blokStarts: comma-separated list of block  
#*							 starts, relative to chromStart
#*
#*    BED is 0-based and half-open. So, "chr 1 100" in a GFF file will
#*    be "chr 0 100" in BED.
#*
#*    Last editted on June 25, 2019
#*
getAnnotationFromBigBedFile <- function(bed_file, has.header=FALSE, 
    sepcial_chrom=NULL) {
    if(!file.exists(bed_file)) stop("Unable to find the file.");
    bed_info <- read.table(bed_file, header=has.header, sep="\t", 
        quote="");
	
    #*  Input file must be in bigBed format (12 columns)
    #*  ---------------------------------------------------------
    if(ncol(bed_info) != 12)
        stop("Bed file must be in bigBed format (12 colums).")
    colnames(bed_info) <- c("chrom", "chromStart", "chromEnd", 
        "name", "score", "strand", "thickStart", "thickEnd", 
        "itemRbg", "blockCount", "blockSizes", "blockStarts");
			
    #*  Gene start position must be less than gene stop position
    #*  since the number increases from 5' to 3' (left to right)
    #*  and gene start cannot be negative
    #*  ---------------------------------------------------------
    gene_start <- as.numeric(bed_info$chromStart);
    gene_stop  <- as.numeric(bed_info$chromEnd);

    if(length(which(gene_start < 0)) > 0 ) 
        stop("Gene start position cannot be less than 0.")
    if(sum(gene_start >= gene_stop) > 0)
        stop("Gene start position greater than stop position.")

    #*  CDS start position must be less than CDS stop position.
    #*  -----------------------------------------------------------
    cds_start <- as.numeric(bed_info$thickStart);
    cds_stop  <- as.numeric(bed_info$thickEnd);

    if(sum(cds_start >= cds_stop) > 0)
        stop("CDS start position greater than stop position.")

	#*  Gene start and stop positions must cover the CDS range.
	#*  Note: the number starts from 5' to 3' (left to right)	
	#*  ---------------------------------------------------------
    if(sum(gene_start > cds_start) > 0 | sum(cds_stop > gene_stop) > 0)
        stop("CDS position is beyond of gene position range.")

    #*  Chromosome names must be in same format
    #*  --------------------------------------------------
    chromosomes <- unique(as.character(bed_info$chrom));    
    total_chroms <- length(chromosomes);

    has.prefix <- grep("chr", chromosomes);
    if(length(has.prefix) > 0 ) {
        if(length(has.prefix) != total_chroms)
            stop("Chromosome name prefix are not uniform.");
        chromosomes <- gsub("chr", "", chromosomes);
        chromosomes <- gsub("^0",  "", chromosomes);
    }

    roman_chroms   <- c(as.character(as.roman(1:100)), "M");
    generic_chroms <- c(as.character(1:100), "X", "Y", "M");
    if(is.null(sepcial_chrom) == FALSE) {
        roman_chroms <- c(roman_chroms, sepcial_chrom);
        generic_chroms <- c(generic_chroms, sepcial_chrom);
    }

    if(sum(chromosomes %in% generic_chroms) != total_chroms &
        sum(chromosomes %in% roman_chroms) != total_chroms)
        stop("Chromosome names are mixed with Roman and digits.")

    #*  Chromosome strand should be in +/- format
    #*  --------------------------------------------------------
    pos_strand <- which(bed_info$strand == "+");
    neg_strand <- which(bed_info$strand == "-");
    if((length(pos_strand) + length(neg_strand)) != nrow(bed_info))
        stop("Chromosome strands should be + or - only.")

    #*  Black size and block start must match to block count
    #*  --------------------------------------------------------

    more_blocks <- which(bed_info$blockCount > 1);
    if(length(more_blocks) > 1) {
        blocks_status <- checkBlocksInBed(bed_info);
        if(blocks_status == FALSE)
            stop("Bad block size or block start definition.");
    }
	
    return (bed_info);
}

#*  =========================================================================
#*
#*    7.  checkBlocksInBed()
#*
#*    Check out the number of block size and block start positions
#*
#*    Argument:
#*
#*        bed_info:	A data frame with 12 columns read from bigBed file.
#*                  It should be checked before call this function and
#*                  must have required standard column names.
#*
#*    Return: Logic, if the number of block size and block start
#*            positions match to block counts.
#*
#*		
#*    Created on June 25, 2019
#*
checkBlocksInBed <- function(bed_info) {

    more_block_lines <- which(bed_info$blockCount > 1);
    if(length(more_block_lines) == 0) {
        message("All rows have one block only.")
        return (TRUE);
    }

    for(a_line in seq_along(nrow(bed_info))) {
        block_number <- as.numeric(bed_info$blockCount[a_line]);
        if(block_number == 1) next;

        sizes_number <- length(strsplit(bed_info$blockSizes[a_line], ",")[[1]]);
        start_number <- length(strsplit(bed_info$blockStarts[a_line],",")[[1]]);
    
        if(block_number != sizes_number || block_number != start_number) {
            message("Bad block definition in line ", a_line, ".");
            return (FALSE);
        }
    }
	
    return (TRUE);
}


#*  =========================================================================
#*    
#*    8.  getSequenceFromMutipleFastaFiles()
#*    
#*    Get genome sequence data from multiple fasta files.
#*   
#*    Arguments:
#*
#*    file_path:  character vector, path to the directory where the fast
#*                files exist 
#*    file_type:  character vector, valid types are "fa", "FA", "fasta", 
#*                and "FASTA"
#*                    
#*    Return:
#*
#*    A data frame containing sequence data where rows are for each 
#*    chromosome and columns are for chromosome names and sequences
#*
#*
#*    Last edited on June 25, 2019
#*
getSequenceFromMutipleFastaFiles <- function(file_path, file_type) {

    if( (file_type %in% c("fa", "FA", "fasta", "FASTA")) == FALSE)
        stop("Incorrect file type defined.")

    fa_files <- list.files(path=file_path, pattern=file_type, full.names=TRUE);
    if(length(fa_files) == 0) stop("No fasta file found.");

    total_chroms <- length(fa_files);
    seq_data <- matrix(rep(" ", 2*total_chroms), ncol=2);
    colnames(seq_data) <- c("chromosome", "sequences")

    for(a_file in 1:length(fa_files)) {
        fa_contents <- readLines(fa_files[a_file], n=-1)
        fa_contents <- gsub("\n", "", fa_contents)
        
        num_chrom <- grep(">", fa_contents);
        if(length(num_chrom) > 1)
            stop("More than one chromosome found in ", fa_files[a_file]);
	
        chromosome <- sub(">", "", fa_contents[1]);
        sequences <- paste(c(fa_contents[2:length(fa_contents)]), collapse="");

        seq_data[a_file, 1] <- chromosome;
        seq_data[a_file, 2] <- sequences; 
    }

    return (as.data.frame(seq_data));
}


#*  =======================================================================
#*
#*    9.  getSequenceFromOneFastaFile()
#*
#*    Get genome sequence data from one fasta files which may contain
#*    sequence for one or more than one chromosome.
#*
#*    Argument:
#*
#*        file_name:  character vector, name (and path) of a fasta file
#*                    
#*    Return:
#*
#*    A data frame wiht row(s) for chromosome(s) and columns for  
#*    chromosome names and sequences.
#*
#*    Last edited on June 25, 2019
#*
getSequenceFromOneFastaFile <- function(file_name) {

    if(!file.exists(file_name)) stop("File does not exist.");

    fa_contents <- readLines(file_name, n=-1);
    fa_contents <- gsub("\n", "", fa_contents);

    chrom_lines <- grep(">", fa_contents);
    if(length(chrom_lines) == 0)
        stop("Unable to find chromosome name line(s).");

    total_chroms <- length(chrom_lines);
    seq_data <- matrix(rep(" ", 2*total_chroms), ncol=2);
    colnames(seq_data) <- c("chromosome", "sequences");
    
    if(length(chrom_lines) > 1) {
        seq_data[,1] <- gsub(">", "", fa_contents[chrom_lines]);
		
        start_lines <- chrom_lines + 1;	
        end_lines <- chrom_lines[2:length(chrom_lines)] - 1;
        end_lines <- c(end_lines, length(fa_contents));
		
        for(a_chrom in 1:length(chrom_lines)) {
            line_number <- start_lines[a_chrom]:end_lines[a_chrom];
            seq_content <- c(fa_contents[line_number]);
            seq_data[a_chrom, 2] <- paste(seq_content, collapse="");    
        }
    } else { 
        seq_data[1, 1] <- fa_contents[1];
        seq_data[1, 2] <- paste(fa_contents[-1], collapse="");
    }

    return (as.data.frame(seq_data));
}

#*  =======================================================================
#*
#*    10.  getUTRInfoDefinedByBed()
#*
#*    Extract genomic sequence for 5'UTR and 3'UTR defined in bed file.  	 
#*    The sequences and their genomic coordinates will be from forward  
#*    (+) strand. In downstream analysis, such as UORF and start codon  
#*    identification, the sequence and coordinates should be modified 
#*    for all features on reverse (-) strand.
#*
#*    Note:
#*  
#*    Genomic sequence is held by character vector and index will start 
#*    from 1. BED is 0-based and half-open. So, "chr 1 100" in GFF file 
#*    will be "chr 0 100" in BED.
#*
#*
#*    Arguments:
#*
#*        DNA_seq:  A data frame with rows for chromosome(s) and 
#*                  columns for chromosome name(s) and sequence.
#*
#*        bed_info:	A data frame with contents same as bigBed file.   
#*                  The data frame should be the one returned by the  
#*                  function call to getAnnotationFromBigBedFile() 
#*                  with column names same as the standard column 
#*                  names of bigBed file.
#*
#*            1).   chrom:      chromosome name of each feature
#*            2).   chromStart: start position of each feature
#*            3).   chromEnd:   end position of each feature
#*            4).   name:       gene name of each feature
#*            5).   score:      used for graphic display only
#*            6).   strand:     chromosome strand of each feature
#*            7).   thickStart: starting position drawn thickly
#*            8).   thickEnd:    ending position drawn thickly
#*            9).   itemRgb:     used for graphic display only
#*            10).  blockCount:  number of blocks of each feature
#*            11).  blockSizes:  comma-separated list of block sizes
#*            12).  blokStarts:  comma-separated list of block  
#*                               starts,relative to chromStart
#*
#*    Return:		
#*
#*        Data frame with 7 columns for:
#*
#*        chromosome:  chromosome name of each UTR 
#*        start_pos:   start position of each UTR 
#*        end_pos:     end position of each UTR 
#*        strand:      strand of each UTR
#*        locus:       gene name of each UTR
#*        sequence:    DNA sequence of each UTR
#*        type:        type of each fragment
#*
#*    Last edited on June 25, 2019
#*
getUTRInfoDefinedByBed <- function(DNA_seq, bed_info) {

    #*  Chromosomes and strands have no change for each fragment
    #*  --------------------------------------------------------
    chromosome  <- as.character(bed_info$chrom);
    strand      <- as.character(bed_info$strand);
    UTR_loc     <- as.character(bed_info$name);

    #*  Each gene will have two UTRs (5_UTR and 3_UTR) and the
    #*  order should be reversed for reverse strand
    #  ------------------------------------------------------
    UTR5_type <- rep("5_UTR", length(chromosome));
    UTR3_type <- rep("3_UTR", length(chromosome));

    rev_strand <- which(strand == "-");
    UTR5_type[rev_strand] <- "3_UTR";
    UTR3_type[rev_strand] <- "5_UTR";

    #  5' side of each feature are temporary set to 5_UTR. Since
    #  the cds are 0-based and half open, the cds start position
    #  will be exactly the last base pair position of UTR and cds
    #  end position plus 1 will be start of 3'UTR
    #  -----------------------------------------------------------
    UTR5_start <- bed_info$chromStart + 1;
    UTR5_end   <- bed_info$thickStart;

    UTR3_start <- bed_info$thickEnd + 1;
    UTR3_end   <- bed_info$chromEnd;

    #*  Extract sequences for each UTR
    #*  ------------------------------------------------------
    seq_chrom <- as.character(DNA_seq[,1]);
    UTR5_seq  <- rep(" ", length(chromosome));
    UTR3_seq  <- rep(" ", length(chromosome));

    for(a_row in 1:nrow(bed_info)) {
        #*  Skip the bed line without feature (CDS).
        #*  --------------------------------------------------
        if(bed_info$blockCount[a_row] == 0) next;

        #*		Chromosome for this bed line
        #*		-----------------------------------------------
        bed_chrom <- as.character(bed_info$chrom[a_row]);
        seq_row   <- which(seq_chrom == bed_chrom);

        UTR5_seq[a_row] <- substr(as.character(DNA_seq[seq_row,2]),
                            UTR5_start[a_row], UTR5_end[a_row])
        UTR3_seq[a_row] <- substr(as.character(DNA_seq[seq_row,2]),
                            UTR3_start[a_row], UTR3_end[a_row])
    }

    UTR_5 <- data.frame(chromosome=chromosome, start_pos=UTR5_start, 
                end_pos=UTR5_end, strand=strand, locus=UTR_loc,  
                sequence=UTR5_seq, type=UTR5_type );

    UTR_3 <- data.frame(chromosome=chromosome, start_pos=UTR3_start, 
                end_pos=UTR3_end, strand=strand, locus=UTR_loc, 
                sequence=UTR3_seq, type=UTR3_type );

    return (rbind(UTR_5, UTR_3));
}

#* =======================================================================
#*
#*    11.  checkChromosomeInfo()
#*
#*    Check if chromosome information in bed file and fasta sequence
#*    match including chromosome names and length.
#*
#*    Arguments:
#*
#*        DNA_seq:	A data frame with rows for chromosome(s) and 
#*                  columns for chromosome name(s) and sequence(s).
#*
#*        bed_info:	A data frame with contents same as bigBed file.   
#*                  The data frame should be the one returned by the  
#*                  function call to getAnnotationFromBigBedFile() 
#*                  with column names same as the standard column 
#*                  names of bigBed file.
#*
#*                1).   chrom:      chromosome name of each feature
#*                2).   chromStart: start position of each feature
#*                3).   chromEnd:   end position of each feature
#*                4).   name:       gene name of each feature
#*                5).   score:      used for graphic display only
#*                6).   strand:     chromosome strand of each feature
#*                7).   thickStart: starting position drawn thickly
#*                8).   thickEnd:   ending position drawn thickly
#*                9).   itemRgb:    used for graphic display only
#*                10).  blockCount: number of blocks of each feature
#*                11).  blockSizes: comma-separated list of block sizes
#*                12).  blokStarts: comma-separated list of block  
#*                      starts,relative to chromStart
#*
#*    Return:  TRUE if chromosome name(s) and length(s) in both bed 
#*             info and fasta info match, otherwise Error message.
#*
checkChromosomeInfo <- function(DNA_seq, bed_info) {

    #*  Input data check. Chromosome names in bedinfo must be
    #*  included in genomicsequences
    #*  -----------------------------------------------------------
    seq_chrom <- unique(as.character(DNA_seq[,1]));
    bed_chrom <- unique(as.character(bed_info$chrom));
    
    if(sum(bed_chrom %in% seq_chrom) != length(bed_chrom))
        stop("Some chromosomes are not found in sequence.")

    #*  Bed range must be in the sequence range
    #*  -----------------------------------------------------
    for(a_chr in 1:length(bed_chrom)) {
        seq_row <- which(DNA_seq$chromosome == bed_chrom[a_chr])
        total_nt <- nchar(as.character(DNA_seq[seq_row,2]));

        chr_row  <- which(bed_info$chrom == bed_chrom[a_chr]);
        last_pos <- max(bed_info$chromEnd[chr_row]);

        if(last_pos > total_nt)
            stop("Bed is beyond chromosome", bed_chrom[a_chr]);
    }
	
    #*	No return if anything is wrong due to the stop().
	     (TRUE);
}


#*  ===================================================================
#*
#*    12.  getUTRSequenceFromFastaFile()
#*
#*    Read UTR sequence from fasta file, e.g., 5UTR.fa downloaded from 
#*    https://www.pombase.org/downloads/utr
#*
#*    A typical record in fast format is much like below:
#*
#*    >SPAC212.06c|SPAC212.06c.1|18558|18974|1||I|protein_coding| 
#*    DNA helicase in rearranged telomeric region, truncated|
#*    CTACACATTACGCTGAGAGGTAAAATACTCTGACAACATTCGTTCGATTGTATAAAACAA
#*    AATCCAGCCGAAACGATTGTTGTCAGTAATCAAGATTACGATCTAAATTGAGTACCAAGA
#*    CAAAACGAAATGGTTAAAAAGTTAAAGTCGTTTTTGTATGGACACAATTTCTATAAAATA
#*    GACATGAGTAAAATCTCGCTATTTGTTTGTTATTGTGGAATAATGAAGAGTCATGGGAGA
#*    TGAATGTTGTAAACGATGGCATAGAATTGGTAACGAAAAGTGAAATCGTTGGGATCAACT
#*    ATTTCAGTATTTTGTTTAAAGAAAATGTTGAACTCGACAAGTAATGAGAGGTGGTGCTTT
#*    CGTTAAATAATGAGTGGTGGTTACGGTTATACAGGATATGATATGTGTATGGTGAGA
#*
#*    The headers are set out like so:
#*
#*    >SPCC338.13|1351474|1354474|-1|cog4|III|protein_coding|Golgi transport 
#*            complex subunit Cog4 (predicted)|
#*     ^          ^       ^       ^  ^    ^   ^              ^
#*     |          |       |       |  |    |   |              \_ Description
#*     |          |       |       |  |    |   \_ Feature type
#*     |          |       |       |  |    \_ Chromosome
#*     |          |       |       |  \_ Gene Name  
#*     |          |       |       \_ Strand
#*     |          |       \_ Gene End
#*     |          \_ Gene Start
#*     \_ Systematic_ID
#*
#*    Note: 
#*
#*    Headers for each feature in example above are in one row 
#*    The number of “|” separator may be more than 1
#*    -------------------------------------------------------------------
#*
#*    Arguments:
#*
#*        fast_file:  character vector, name (and path) of fasta file
#*        headers:    character vector, definition of each header field.
#*
#*    Return:
#*
#*        A data frame with 2 columns for chromosome names and sequence
#*        of each UTR
#*
getUTRInfoFromFastaFile <- function(fasta_file, header_name) {

    # Read all lines first.
    # -------------------------------------------------------
    fasta_data <- readLines(fasta_file, n=-1)

    head_lines <- grep(">", fasta_data);
    header_content <- fasta_data[head_lines];

    ##  Extract header information
    ##  -------------------------------------------------------
    header_num <- length(header_name);
    header_items <- matrix(rep(" ", header_num*length(head_lines)), 
                       ncol=header_num)
    for(index in 1:length(header_content))
    {
        headers <- unlist(strsplit(header_content[index], split="\\|"));
        headers <- headers[-which(headers == "")];
        headers[1] <- sub(">", "", headers[1])

        if(length(headers) != header_num) 
            stop("Header lines are not uniform.");

        header_items[index,] <- headers;
    }
    colnames(header_items) <- header_name;
 

    #    Concatenate sequence lines for each UTR
    #    --------------------------------------------------------
    start_lines <- head_lines + 1;
    end_lines <- head_lines[2:length(head_lines)] - 1;
    end_lines <- c(end_lines, length(fasta_data));

    sequences <- rep(" ", length(head_lines))
    for(index in 1:length(head_lines))
    {
        utr <- c(fasta_data[start_lines[index]:end_lines[index]])
        utr <- gsub("\n", "", utr);
        sequences[index] <- paste(utr, collapse="");
    }

    UTR <- cbind(header_items, sequences);

    return (as.data.frame(UTR));
}

#*  =======================================================================
#*
#*    13.  covertToRNASequence()
#*
#*    Convert a fragment of DNA sequence to mRNA sequence, For fragment
#*    on forward strand, simply convert "T" to "U". For fragment on  
#*    reverse strand first get the reverse strand sequence, then reverse 
#*    to 5' -> 3' order and replace "T" with "U"
#*
#*    Arguments:
#*
#*        DNA_seq:  character vector, a fragment of DNA sequences
#*        strand:   character, either "+" or "-"
#*
#*    Return:
#*
#*        A character vector containing series of "A", "U", "G", "C".
#*
#*    Last edited on June 25, 2019
#*
covertToRNASequence <- function(DNA_seq, strand) {

    if(strand != "+" & strand != "-")
        stop("strand must be + or -.");

    DNA_seq <- toupper(DNA_seq);
    if(is.DNA.sequence(DNA_seq) == FALSE)
        stop("Sequences has bases other than A, T, G, C");

    if(strand == "-") {
        DNA_seq <- getComplementarySequence(DNA_seq);
        DNA_seq <- getReversedSequence(DNA_seq);
    }

    mRNA_seq <- gsub("T", "U", DNA_seq);

    return (mRNA_seq);
}


#*  =======================================================================
#*
#*    14.	getReversedSequence()
#*		
#*    Reverse a fragment of sequences withput loading other packages.
#*    
#*    Argument:
#*
#*    seq_fragment:  character vector, a fragment of DNA sequence.
#*
#*    Return:
#*
#*        Character vector, same base pairs as input DNA sequences
#*        in reversed order (ATGC -> CGTA)      
#*
getReversedSequence <- function(seq_fragment) {

    bases <- strsplit(seq_fragment, "")[[1]];
    bases <- bases[length(bases):1];
    bases <- paste(bases, collapse="");

    return (bases);
}

#*  ===================================================================
#*
#*    15.	getComplementarySequence()
#*
#*    Get complementary sequence for a fragment of DNA sequence
#*
#*    Argument:
#*
#*        seq_fragment:  character vector, a fragment of DNA sequence.
#*
#*    Return:
#*
#*        character vector, the complementary sequence of the argument 
#*        where A -> T, G -> C, T -> A, and C -> G.
#*
#*    Last edited on June 25, 2019
#*
getComplementarySequence <- function(seq_fragment) {

    cSequence <- toupper(seq_fragment);

    cSequence <- gsub("A", "1", cSequence);
    cSequence <- gsub("T", "2", cSequence);
    cSequence <- gsub("G", "3", cSequence);
    cSequence <- gsub("C", "4", cSequence);

    cSequence <- gsub("1", "T", cSequence);
    cSequence <- gsub("2", "A", cSequence);
    cSequence <- gsub("3", "C", cSequence);
    cSequence <- gsub("4", "G", cSequence);

    return (cSequence);
}

#* ===================================================================
#*
#*    16.  is.DNA.sequence()
#*
#*    Check if a DNA sequence fragment contains any bases other than 
#*    "A", "T", "G", or "C"
#*
#*    Argument:
#*
#*        DNA_seq:  character vector, a fragment of DNA sequence.
#*
#*    Return:
#*
#*        Logic, TRUE if argument is DNA sequence, otherwise FALSE
#*
#*    Last edited on June 25, 2019
#*

is.DNA.sequence <- function(DNA_seq) {
    if(is.character(DNA_seq) == FALSE) {
        message("Sequences must be character vector.")
        return (FALSE);
    }
    sequences <- toupper(DNA_seq);

    DNA_bases <- c("A", "T", "G", "C");
    seq_bases <- strsplit(DNA_seq, split="")[[1]];

    if(sum(seq_bases %in% DNA_bases) != length(seq_bases)) {
        message("DNA sequence bases must be A, T, G or C.")
        return (FALSE);
    } else {
        return (TRUE);
    }
}


#*  ===================================================================
#*
#*    17.  getCDSSequenceDefinedByBed()
#*
#*    Extract genomice sequence for CDS defined in bed file. The CDS 
#*    include all length from first thickStart to last thickEnd (e.g.,
#*    exclude 5'UTR and 3'UTR regions)
#*
#*    Arguments:
#*
#*        DNA_seq:   A data frame with rows for chromosome(s) and  
#*                   columns for chromosome name(s) and their 
#*                   sequence(s).
#*
#*        bed_info:  A data frame with contents same as bigBed file. The  
#*                   data frame should be the one returned by function 
#*                   call to getAnnotationFromBigBedFile() and column  
#*                   names are same as the standard column names of 
#*                   bigBed file format.
#*
#*            1).   chrom:       chromosome name of each feature
#*            2).   chromStart:  start position of each feature
#*            3).   chromEnd:    end position of each feature
#*            4).   name:        gene name of each feature
#*            5).   score:       used for graphic display only
#*            6).   strand:      chromosome strand of each feature
#*            7).   thickStart:  starting position drawn thickly
#*            8).   thickEnd:    ending position drawn thickly
#*            9).   itemRgb:     used for graphic display only
#*            10).  blockCount:  number of blocks of each feature
#*            11).  blockSizes:  comma-separated list of block sizes
#*            12).  blokStarts:  comma-separated list of block  
#*                               starts,relative to chromStart
#*
#*    Return:  Data frame with 7 columns for:
#*
#*        chromosome:  chromosome name of each UTR 
#*        start_pos:   start position of each UTR 
#*        end_pos:     end position of each UTR 
#*        strand:      strand of each UTR
#*        locus:       gene name of each UTR
#*        sequence:    DNA sequence of each UTR
#*        type:        type of each fragment
#*
#*		Last edited on June 25, 2019
#*
getCDSInfoDefinedByBed <- function(DNA_seq, bed_info) {

    #*  Chromosomes and strands have no change for each fragment
    #*  --------------------------------------------------------
    CDS_chrom  <- as.character(bed_info$chrom);
    CDS_strand <- as.character(bed_info$strand);
    CDS_locus  <- as.character(bed_info$name);
    CDS_type   <- rep("CDS", length(CDS_chrom));
 
    CDS_start <- bed_info$thickStart + 1;
    CDS_end   <- bed_info$thickEnd;

    #*  Extract sequences for each CDS
    #*  -----------------------------------------------------------
    seq_chrom <- as.character(DNA_seq[,1]);
    CDS_seq   <- rep(" ", length(CDS_chrom));
    bed_chrom <- as.character(bed_info$chrom);

    for(a_row in 1:nrow(bed_info)) {
        seq_row   <- which(seq_chrom == bed_chrom[a_row]);
        CDS_seq[a_row] <- substr(as.character(DNA_seq[seq_row,2]),
            CDS_start[a_row], CDS_end[a_row]);
    }

    CDS <- data.frame(chromosome=CDS_chrom, start_pos=CDS_start, 
                end_pos=CDS_end, strand=CDS_strand, locus=CDS_locus, 
                sequence=CDS_seq, type=CDS_type );

    return (CDS);
}

#* =======================================================================
#*
#*    18.	writeORFInfoToFile()
#*
#*    The main function to screen ORF information. 
#*
#*    For each UTR/CDS, the DNA sequence is first converted to mRNA 
#*    sequence then splited to codon list starting from the 1st, 2nd, 
#*    and 3rd base (frames).
#*
#*    After having the codon list, start and stop codons are searched 
#*    (each UTR may have zero or more than one ORF) and the start and 
#*    stop positions of ORF, length of ORF, distance to UTR start and 
#*    end positions are calculated.
#*
#*    Information of all ORF are written to a tab-delimited text file.
#*
#*    Arguments: 
#*
#*        seq_info:      Data frame, UTR/CDS info including chromosome,  
#*                       start and stop position, strand, sequence, and
#*                       type (CDS, 5_UTR or 3_UTR).
#*        start_codons:  Character vector, start codons.
#*        stop_codons:   Character vector, stop codons.
#*        out_file:      Character vector, name (and path) of output
#*                       file to write ORF information.
#*
#*    Return:   None. Results are written to file.
#*
#*    Last edited on June 27, 2019
#*
writeORFInfoToFile <- function(seq_info, start_codons, stop_codons, 
    out_file) {
    #*  Write column headers to file first
    #*  -------------------------------------------------
    headers <- c("gene_name", "chromosome", "seq_start",
            "seq_stop", "strand", "orf_id", "orf_start", 
            "orf_stop", "start_codon", "context",   
            "orf_length", "dis_from_cap", "dis_to_MUG");
    write(headers, file=out_file, ncolumns=length(headers), 
        sep="\t", append=FALSE);
	
    for(a_seq in 1:nrow(seq_info))  {
        message(a_seq, ": ", seq_info$locus[a_seq]);
		
        #*  Convert the DNA sequence in seq_info to RNA sequence 
        #*  -----------------------------------------------------
        the_seq  <- as.character(seq_info$sequence[a_seq]);
        strand   <- as.character(seq_info$strand[a_seq]);
        mRNA_seq <- covertToRNASequence(the_seq, strand);
        mRNA_length <- nchar(mRNA_seq);

        #*  Get ORF information from a sequence fragment. Repeat 
        #*  reading codons by starting at the 1st, 2nd, and 3rd
        #*  base from the 5' of sequence. Sequence with less than 
        #*  3 nucleotides will be skipped.
        #*  -----------------------------------------------------
        for(start_at in 1:3) {
            if((mRNA_length - start_at) < 2) next;
            codon_list <- getCodonFromSequence(mRNA_seq, start_at);

            start_index <- getCodonIndex(codon_list, start_codons);
            if(length(start_index) == 0 ) next;

            stop_index  <- getCodonIndex(codon_list, stop_codons)
            if(length(stop_index) == 0) stop_index <- 0;

            codon_pair <- setCodonIndexPair(start_index, stop_index);
            orf <- screenORF(seq_info[a_seq,], codon_list, codon_pair,  
                        start_at, mRNA_seq);
			
            write.table(orf, file=out_file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE, append=TRUE)
        }    
    }
}


#*  =======================================================================
#*
#*    19.  getCodonIndex()
#*
#*    Match codon list to target codons.
#*
#*    Arguments:
#*
#*        codon_list:    Character vector, a serial codons from a
#*                       sequences fragment
#*                
#*        target_codon:  Character vector, one or more start or 
#*                       stop codon.
#*
#*    Returne:
#*
#*        Numeric vector, which target codon a codon in codon list
#*        matches. If target codon matched, zero will be returned.
#*
#*    Last edited on June 26, 2019
#*
getCodonIndex <- function(codon_list, target_codon) {
    validateCodons(codon_list);
    validateCodons(target_codon);

    index <- which(codon_list %in% target_codon );

    ##  This does not need, just want more safe
    index <- sort(index);
  
    return (index);
}


#* =======================================================================
#*
#*    20.  setCodonIndexPair()
#*
#*    Pair start and stop codons in an ORF ranges. For each start codon,
#*    only the first stop codon in its downstream can be used. If there  
#*    is no stop codon found, the end of sequence is used (index 0).
#*      
#*    Argument:
#*
#*        start_codons:  integer vector, start codon index, always 
#*                       greater than 0. Otherwise this function 
#*                       will not be called.
#*        stop_codons:   integer vector, stop codon index,  either 0 
#*                       or integer(s) greater than 0.
#*
#*    Return: 
#*
#*        A matrix containing paired codon index which representing 
#*        ORF range and total number of codons.
#*
#*    Last edited on June 26, 2019
#*
setCodonIndexPair <- function(start_index, stop_index) {
    #* Initialize a matrix for start and stop codon pairs
    #* --------------------------------------------------
    codon_pairs <- cbind(start_index, rep(0, length(start_index)));
    
    #  Stop codon must be in down stream of start codon 
    #  and if there is no stop codon, set index to 0
    #  -------------------------------------------------
    
    for(a_row in 1:length(start_index)) {
        big_end <- which(stop_index > codon_pairs[a_row, 1]);
        if(length(big_end) > 0) 
            codon_pairs[a_row, 2] <- stop_index[big_end[1]];
    }

    return (codon_pairs);
}


#* =======================================================================
#*
#*    21.	screenORF()
#*
#*    Screen a sequence fragment to get start codon, context from -3  
#*    to +4 base positions of start codon, length of ORF, distance from 
#*    the cap, and distance to the uAUG in a list.
#*       
#*    Arguments:
#*
#*        seq_info:    One row of data frame, sequence info including of 
#*                     chromosome, start_pos, end_pos, strand, locus, 
#*                     sequence, and type.
#*        codon_list:  Codons converted from UTR sequence
#*        codon_pair:  Matrix of column 2 for start and stop index  
#*        start_at:    Integer in 1, 2, 3, base where start to read
#*
#*    Return: 
#*
#*        A data frame containing all parameters above
#*
#*    Last edited on June 26, 2019
#*
screenORF <- function(seq_info, codon_list, codon_pair, 
    start_at, mRNA_seq) {
    total_orf <- nrow(codon_pair);
    ORF_info <- data.frame();

    for(a_row in 1:total_orf) {
        #*  Context of ORF
        ##  ----------------------------------------------------
        start_index <- as.numeric(codon_pair[a_row, 1]);
        stop_index  <- as.numeric(codon_pair[a_row, 2]);
        orf_context <- getStartCodonContext(codon_list, start_index, 
                            mRNA_seq, start_at);

        #*		Length of ORF. Stop codon will not be included
        #*		----------------------------------------------------
        last_codon <- stop_index;
        if(last_codon == 0) {  
            last_codon <- length(codon_list);
        } else { last_codon <- stop_index - 1; }
        orf_length  <- length(start_index:last_codon) * 3;

        #*  Positions of ORF
        #*  ---------------------------------------------------
        orf_positions <- getORFPositions(seq_info, start_index, 
            stop_index, start_at) 
        orf_start <- orf_positions[1]; 
        orf_stop  <- orf_positions[2];
        from_cap  <- orf_positions[3];
        to_MUG    <- orf_positions[4];

        #*  ORF name (gene_name, UTR/CDS, num_ORF, dis_from_cap)
        #*  -----------------------------------------------------
        orf_name <- paste(seq_info$locus, seq_info$type, from_cap, 
            a_row, sep="_");
    
        #*  Put every item into data frame
        #*  -----------------------------------------------------
        a_orf <- data.frame(
            gene_name=as.character(seq_info$locus), 
            chromosome=as.character(seq_info$chromosome), 
            seq_start=as.numeric(seq_info$start_pos),
            seq_stop=as.numeric(seq_info$end_pos),
            strand=as.character(seq_info$strand), 
            orf_id=orf_name, 
            orf_start=orf_start, orf_stop=orf_stop,
            start_codon=codon_list[start_index], 
            context=orf_context,   
            orf_length=orf_length, 
            dis_from_cap=from_cap, 
            dis_to_MUG=to_MUG);
       
        ORF_info <- rbind(ORF_info, a_orf);
    }

    return (ORF_info);
}


#*  =======================================================================
#*
#*    22.	getORFPositions()
#*
#*    Calculate genomic positions for a ORF (start, stop, distance
#*    to cap (UTR start), and distance to main AUG codon
#*
#*    Arguments:
#*
#*        seq_info:     information of UTR (one row from data frame) 
#*        start_index:  integer, start codon index in an ORF  
#*        stop_index:   integer, stop codon index in an ORF
#*        start_at:     integer, start point to read codons
#*
#*    Return:
#*
#*        Numeric vector (all integers) for start and stop positions,
#*        distance from cap and to main AUG.
#*
#*    Last edired on June 26, 2019
#*
getORFPositions <- function(seq_info, start_index, stop_index, start_at, 
                        include_stop_codon) {
    #*  Positions for UTR are from UTR in 5' to 3' direction
    #*  -----------------------------------------------------
    seq_start  <- as.numeric(seq_info$start_pos)
    seq_stop   <- as.numeric(seq_info$end_pos);
    seq_length <- length(seq_start:seq_stop);
    seq_strand <- as.character(seq_info$strand);
  
    orf_start <- -1;
    orf_stop <- -1;
    dis_from_cap <- -1;
    dis_to_mAUG <- -1;
    
    #*  Relative positions of ORF in 5' to 3' direction. The
    #*  start position is the first base of start codon and
    #*  stop position is the last base of last coding codon
    #*  ----------------------------------------------------
    start_position <- (start_index - 1) * 3 + start_at;    
    dis_from_cap   <- start_position - 1;
    dis_to_mAUG    <- seq_length - start_position;

    if(stop_index == 0) {
        stop_position <- seq_length;
    } else {
        stop_position <- (stop_index - 1) * 3 + start_at - 1;
    }

    #*  Convert to genomics positions inside of sequence
    #*  -------------------------------------------------
    
    #*  Forward strand, 5' to 3' 
    #*  -------------------------------------------------
    if(seq_strand == "+") {
        orf_start <- seq_start + dis_from_cap;
        orf_stop  <- seq_start + stop_position - 1; 

    #*  Reverse strand, two positions need be switched
    #*  ------------------------------------------------
    } else if (seq_strand == "-") {
        orf_start <- seq_stop - stop_position + 1;
        orf_stop  <- seq_stop - dis_from_cap; 

    #*  Error message, should be never happen
    #*  -------------------------------------
    } else {
        message("Looks like there is undefined strand.");
    }

    orf_positions <- c(orf_start, orf_stop, dis_from_cap, dis_to_mAUG);

    return (orf_positions)
}


#*  =======================================================================
#*
#*    23.  getStartCodonContext()
#*
#*    Build the context (a short sequence fragment from the third 
#*    nucleotide before start codon and the one nucleotide next to 
#*    start codon). There will be total of 7 nucleotides, e.g., 
#*    GACAUGG, AUCAUGC) and the bases are extracted from sequence.
#*    dis_from_cap is used here since it is relative to sequence 
#*    start position.
#*
#*    Arguments:
#*
#*        codon_list:   Character vector, codons from sequence
#*        start_index:  Integer, index of start codon. Always 
#*                      greater than 0.
#*        mRNA_seq:     Character vector, mRNA sequence converted 
#*                      from UTR/CDS sequence
#*        start_at:     Integer, either 1, 2 or 3.
#*
#*    Return:  
#*	
#*        character vector, nucleotide around start codon.  
#*	
#*    Last edited on June 26, 2019
#*
getStartCodonContext <- function(codon_list, start_index, 
    mRNA_seq, start_at) {

    first_3 <- substr(mRNA_seq, 1, 3);
    last_3  <- substr(mRNA_seq, nchar(mRNA_seq) - 2, nchar(mRNA_seq));

    start_codon <- codon_list[start_index];

    #*  the three bases before start codon, if any
    #*  ------------------------------------------
    if(start_index > 1) { 
        prev_bases <- codon_list[start_index-1];
    } else {
        if(start_at == 1) {
            prev_bases <- "";
        } else {
            prev_bases <- substr(mRNA_seq, 1, start_at-1);
        }      
    }

    #*  the last base after stop codon, if any
    #*  -----------------------------------------------
    if( start_index < length(codon_list)) {
        next_base <- substr(codon_list[start_index+1], 1, 1)
    } else {
        last_base <- length(codon_list)*3 + start_at - 1;
        if(last_base == nchar(mRNA_seq) ) {
            next_base <- "";
        } else {
            next_base <- substr(mRNA_seq, last_base+1, 1);
        }
    } 

    orf_context <- paste0(prev_bases, start_codon, next_base);  

    return (orf_context);
}


#* =======================================================================
#*
#*    24. converToGenePredFormat()
#*
#*    Covert ORF table to GenePrep format table. All positions will be
#*    1-based. ORF with 0 codons between start and stop codons will be 
#*    removed.
#*		
#*    Argument:  
#*
#*        ORF_table:    Data frame returned converted from object
#*                      returned by getORFInfo()
#*        is.one.based:	Logic, if position is 1-based. Default TRUE.
#*
#*    Return: 
#*
#*        Data frame with customized columns of GenePrep format.
#*
#*    Last edited on June 26, 2019
#*
converToGenePredFormat <- function(ORF_file, is.one.based=TRUE, 
    min_codon=0) {

    ORF_table <- read.table(ORF_file, header=TRUE, sep="\t",
        quote="", stringsAsFactors=FALSE);

    orf_start <- ORF_table$orf_start;
    if(is.one.based == FALSE) {
        orf_start <- orf_start - 1;
    }

    #*  Per James method, not for others
    #*  ------------------------------------------------------
    gpTable <- data.frame(gene_name=ORF_table$gene_name,
                  orf_ID=ORF_table$orf_id,
                  chromosomes=ORF_table$chromosome,
                  strand=ORF_table$strand,
                  tx_start=orf_start,
                  ts_end=ORF_table$orf_stop,
                  cds_start=orf_start,
                  cds_end=ORF_table$orf_stop,
                  num_exon=rep(1, nrow(ORF_table)),
                  exon_start=paste0(orf_start, ","),
                  exon_end=paste0(ORF_table$orf_stop, ",")
    )

    #*  keep ORFs with at least two codons (one start and one 
    #*  non-stop codon), meaning total of 6 nucleotides and 
    #*  stop_position - start_position must be greater than 5
    #*  (6 - 1 = 5)
    #  -----------------------------------------------------
    if(min_codon > 0){
        rows <- which((ORF_table$orf_stop - ORF_table$orf_start) >= 5);
        if(length(rows) > 0) gpTable <- gpTable[rows,];
    }
 
    return (gpTable);
}


#*
#*  End of RiboPRoR_ORF_Utility.R
#*  ===================================================================