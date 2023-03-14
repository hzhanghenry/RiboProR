#*
#*    File:	RiboProTools_Fp_Density_Center.R
#*
#*    Find the center of footprint density for any transcript to
#*    measure if there are more transcription on 5' uORF region.
#*
#*    Functions implemented:
#*
#*    1.  getFootprintDensityCenter()
#*    2.  getOriginalFPCenter()
#*    3.  getFilteredFPCenter
#*    4.  getTranscriptLength()
#*    5.  getCDSGRanges()
#*
#*    Last revised on August 10, 2018
#*  ________________________________________________________________________
#*  <Fp density><Fp density><Fp density><Fp density><Fp density><Fp density>



#*  =================================================================
#*
#*    1.  The main function for footprint center weighting
#*
#*    Pipeline:
#*
#*    a_site <- FpFraming(bam_file, bed_file, save.file=FALSE);
#*    YeastBedGRange <- import(bedfile, format="BED");
#*    transcript <- transcriptGRanges(YeastBedGRange[1]);
#*    asiteProfile <- getASiteProfile(asites, bamfile, transcript);
#*    fp_center <- getFPCenter(asiteProfile);
#*
#*    Argument::
#*
#*        a_site:       Data frame with 1 column for a-site. Row
#*                      names are read length.
#*        transcripts:  GRanges list of all transcripts.
#*        bam_file:     Character vector, a bam file name (and path).
#*
#*    Return:
#*
#*        A data frame with one column of positive integers for density
#*        center of all transcripts.
#*

getFootprintDensityCenter <- function(a_sites, transcripts, bam_file, 
		weight=FALSE) {
    centers_list <- mclapply(1:length(transcripts), 
        function(i) {
            asiteProfile <- getASiteProfile(a_sites, bam_file, transcripts[i]);
            trx_length <- sum(width(transcripts[i]$blocks[[1]]));
 
            if(weight == FALSE) {
                fp_center <- getOriginalFPCenter(asiteProfile, trx_length);
            } else {
                fp_center <- getFilteredFPCenter(asiteProfile);
            }
            return (fp_center)
        }
    )
    densisy_centers <- unlist(centers_list);
	
    return (data.frame(densisy_centers, row.names=transcripts$name));

}


#*  =========================================================================
#*
#*    2. getOrignalFPCenter()
#*
#*    Find center of footprints for a transcript without filtering
#*
#*    Arguments:
#*
#*        asiteProfile:  Vector of positive integer, counts of a-site
#*                       on each nucleotide position of a transcript
#*        trx_length:    Positive integer, length of the transcript
#*
#*    Return:
#*
#*        Positive float number, ratio of length of the 5' end of transcript
#*        that has half of total counts divided by transcript length
#*

getOriginalFPCenter <- function(asiteProfile, trx_length) {

    the_center <- NA;
    total <- sum(asiteProfile);

    if(total > 0) {
        half_of_total <- total/2;
        count <- 0;
        for(a_site in seq_along(asiteProfile)) {
            count <- count + asiteProfile[a_site];
            if(count >= half_of_total) break;
        }
        the_center <- a_site/trx_length;
    }
	
    return (the_center)
}


#*  =================================================================
#*
#*    3.  getFilteredFPCenter()
#*
#*    Find center of footprints for a transcript with zero count position
#*    being filtered out.
#*
#*    Arguments:
#*
#*        asiteProfile:	Vector of positive integer, counts of a-site
#*                      on each nucleotide position of a transcript
#*
#*    Return:
#*
#*        Positive float number, length of the 5' end of transcript
#*        that has half of total counts divided by transcript length
#*

getFilteredFPCenter <- function(asiteProfile) {

    zeros <- which(asiteProfile == 0)
    asiteProfile <- asiteProfile[-zeros];
	
    total <- sum(asiteProfile);
    the_center <- NA;
    if(total > 0) {
        half_of_total <- total/2;
        count <- 0;
        for(a_site in seq_along(asiteProfile)) {
            count <- count + asiteProfile[a_site];
            if(count >= half_of_total) break;
        }
        the_center <- a_site/length(asiteProfile);
    }
	
    return (the_center)
}


#*  =================================================================
#*
#*    4.  getTranscriptLength()
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
    return (data.frame(counts, row.names=transcripts$name));
}


#*  =================================================================
#*
#*    5.  getCDSGranges()
#*
#*    Generate GRanges object for cds regions from transcripts 
#*    defined with bed file
#*
#*    Argument:
#*
#*        bed_file:  character vector, name (and path) of a bed file
#*                   defining transcripts
#*
#*    Return:  GRanges object for cds region only 		
#*
getCDSGRanges <- function(bed_file) {

    file_ext <- toupper(sub(".{1,}\\.", "", bed_file));
    if(file_ext != "BED") stop("Need bed file");
	
    trx <- import(bed_file, "BED");
    cds <- trx
    start(cds) <- start(trx$thick)
    end(cds) <- end(trx$thick)

    cds$blocks <- IRanges::IRangesList(mclapply(1:length(cds), 
        function(i) {
            cds_block <- IRanges::IRanges(start=1, 
                end=end(cds[i]) - start(cds[i]) + 1);
            return (cds_block);
        })
    );
	
    return (cds);
}

##  End of RiboProTools_Fp_Density_venter.R
##  ===========================================================================
