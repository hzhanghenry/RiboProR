#*
#*    File:  RiboProfiling_Plot_Util_Source.RGB
#*
#*    Plot Utilities for Presenting Riboseq Profiling Data
#*
#*    Function implemented:
#*
#*    1)   getCountsData()
#*    2)   fpCountsBoxPlot()
#*    3)   plotReplicates()
#*    4)   plotTranslationEfficiency()
#*    5)   getPlotColors()
#*    6)   getDefaultColors()
#*    7)   fpFramePlot()
#*    8)   plotMetageneFrames
#*    9)   getRegionFrames()
#*    10)  getPrettyLabels()
#*    11)  plotMetageneCounts() 		
#*    12)  plotHeatmap()
#*    13)  plotCorrelationHeatmap()		
#*    14)  plotRedBlueCorrelationImage()
#*    15)  plotCorrelationMatrix()
#*    16)  plotMultiMetageneCounts()
#*    17)  getMetageCountsMatrix()
#*
#*  Last revised on January 08, 2019
#*  __________________________________________________________________
#*  <Plot_Util><Plot_Util><Plot_Util><Plot_Util><Plot_Util><Plot_Util>



#*  ==================================================================
#*
#*    1)  getCountsData(dds)
#*
#*    Get raw counts data or normalized counts from DESeqDataset object 
#*    for visualization.
#*
#*    Arguments:
#*
#*        dds:        A DESeqDataSet object which holds the raw counts. 
#*        normalize:  Character vector, methods for normalization 
#*                    including of "size", "fpm", and "fpkm". 
#*                    If normalizing for fpkm, basepairs column must
#*                    exist in mcols of dds. Could be set as NULL
#*                    for no normalization. (default: "size") 
#*
#*    Return:  A matrix with normalized count data
#*
getCountsData <- function(dds, normalize=c("size", "fpm", "fpkm")){

    if(is.null(normalize)) {
        count_data <- counts(dds, normalize=FALSE);
    } else {
        normalize <- tolower(normalize[1]);
		
        if(normalize == "size") {
            if(is.null(sizeFactors(dds)))
                dds <- estimateSizeFactors(dds);
            count_data <- DESeq2::counts(dds, normalize=TRUE);
        } else if (normalize == "fpm") {
            count_data <- fpm(dds);
        } else if (normalize == "fpkm") {
            if(length(grep("basepairs", colnames(mcols(dds)))) == 0 && 
                length(grep("avgTxLength", colnames(mcols(dds)))) == 0)
                stop("Required columns for fpkm is not found.")
            count_data <- fpkm(dds);
        } else { 
            stop(paste("Normalization method", normalize, 
                "not supported")); 
        }
    }
	
    return (count_data);
}


#*  ==================================================================
#*
#*    2)  fpCountsBoxPlot()
#*
#*    Box plot for basic description of foot print (fp) counts data
#*
#*    Arguments:
#*
#*        dds:         A DESeqDataSet object with size factor calculated   
#*                     or basepairs column presented in mcols.	
#*        normalize:   Character vector, normalization methods including 
#*                     of "size", "fpm", and "fpkm". Could be set to null 
#*                     for no normalization.
#*        label_area:  Positive integer, height of area at botom of plot  
#*                     for x axis labels.
#*        main_text:   Character vector, text for title of the plot.
#*
#*    Return:  None. Plot figure only
#*
fpCountsBoxPlot <- function(dds, normalize=c("size", "fpm", "fpkm"),
    label_area=2, main_text="Distribution of fp Counts") {

    count_data <- getCountsData(dds, normalize);
    count_data <- log2(count_data + 1);

    plot_colors <- getPlotColors(ncol(count_data));

    par(mai=c(label_area,1,1,1))
    boxplot(count_data, las=2, notch=TRUE, col=plot_colors);
    title(main_text);
}


#*  ==================================================================
#*
#*    3)  plotReplicates()
#*
#*    Scattern plot of two replicate samples in dataset.
#*
#*    Argument(s):
#*
#*        dds:          A DESeqDataSet object with size factor calculated  
#*                      or basepairs column presented in mcols.
#*        normalize:    Character vector, normalization methods including 
#*                      of "size", "fpm", and "fpkm".
#*        replicates:   Positive integer vector of length 2. The columns for 
#*                      replicates to be plotted. 
#*        show.cor:     Logic, if plot corelation coefficiency and the 
#*                      regression line.
#*        is.log2:      Logic, is the data log2 transformed
#*        point_color:  Character vector of R colors names, color for the
#*                      point background. 
#*        line_color:   Character vector, colors for point or regression line.
#*        x_label:      Character vector, text for x axis label
#*        y_label:      Character vector, text for y axis label
#*        main_text:    Character vector, title of the plot
#*
#*    Return:  None. Plot figure only
#*	

plotReplicates <- function(dds, normalize=c("size", "fpm", "fpkm"),
    replicates=c(1,2), show.cor=TRUE, is.log2=FALSE,
    point_color="grey", line_color="red",
    x_label="replicate 1", y_label="replicate 2",
    main_text="Distribution of fp Counts") {

    count_data <- getCountsData(dds, normalize);
    count_data <- log2(count_data + 1);
	
    x <- as.numeric(count_data[, replicates[1]]);
    y <- as.numeric(count_data[, replicates[2]]);
    plot(x, y, pch=21, bg=point_color, xlab=x_label, ylab=y_label);
    abline(lm(x ~ y), col=line_color, lwd=2);
    title(main_text, cex=0.8);
	
    if(show.cor == TRUE) {
        cor_res <- cor.test(x, y);
        r <- signif(cor_res$estimate, digits=4);
		
        text_x <- max(x)*0.75;
        text_y <- max(y)*0.25;
        text(text_x, text_y, paste("r =", r));
    }
}


#*  =======================================================================
#*
#*    4)  plotTranslationEfficiency()
#*
#*    Scatterplot with too translational efficiency (TE) values, e.g,
#*    TE of wildtype and TE of mutant. Also mark the TE changes higher
#*    than ratio_level as red and the TE Change below ratio_level as 
#*    blue, and draw a regression line. 
#*
#*    Arguments:
#*
#*        dds:          A DESeqDataSet object returned from funtion of  
#*                      extractEfficiencyChange().
#*        ratio_level:  Positive numeric, threshold to change point colors.
#*        x_label:      Character vector, text for x axis label.
#*        y_label:      Character vector, text for y axis label.
#*        title_text:   Character vector, text for title.
#*        x_pos:        Integer, x coordinate for text plot of pearson's r
#*        y_pos:        Integer, y coordinate for text plot of pearson's r
#*
#*    Return:  None. Plot figure only
#*
plotTranslationEfficiency <- function(dds, ratio_level=1, 
    x_label="TE of Wild Type", y_label="TE of Mutant", 
    title_text="Translation Efficiency Of Mutant",
    x_pos=0, y_pos=12) {

    #*	dds object must be from TE analysis
    #*	-----------------------------------
    if(resultsNames(dds)[3] != "condition_Ribo_vs_mRNA")
        stop("The dds object must be from TE analysis.")
	
    control_TE <- DESeq2::results(dds, contrast=c("condition","Ribo","mRNA"));
    mutant_TE <- DESeq2::results(dds, 
        list(c(resultsNames(dds)[3], resultsNames(dds)[4])));

    mut_te <- as.numeric(mutant_TE[,2]);
    wt_te  <- as.numeric(control_TE[,2]); 
    TE_fc  <- mut_te - wt_te;
	
    colors <- rep("black", length(mut_te));
    colors[TE_fc >= ratio_level] <- "red";
    colors[TE_fc <= -1*ratio_level] <- "blue";
    r <- cor(mut_te, wt_te); 
	
    par(mfrow=c(1,1), mai=c(1,1,1,1));
    plot(wt_te, mut_te, xlab=x_label, ylab=y_label, 
        pch=21, bg="grey", col=colors);
    abline(lm(mut_te ~ wt_te), col="darkgreen", lwd=2);
    text(x_pos, y_pos, paste("r =", signif(r, digits=4)));
    title(title_text, cex.main=0.8);
}


#*  =======================================================================
#*
#*    5)  getPlotColors()
#*
#*    Get plot colors based on number of samples (columns) in plot data.
#*    Maximum number of colors is 24.
#*
#*    Argument: 
#*
#*        num_colors: positive integer, total number of colors.
#*
#*    Return:  Vector of R colors with length of num_colors.
#*
getPlotColors <- function(num_colors) {

    plot_colors <- getDefaultColors();
	
    if(num_colors > length(plot_colors)) {
        plot_colors <- rainbow(num_colors);
    } else {
        plot_colors <- plot_colors[1:num_colors];
    }
    return (plot_colors);
}


#* =======================================================================
#*
#*    6)  getDefaultColors()
#*
#*    Get default plot colors (total of 24) 
#*	
#*    Argument: None
#*    Return:   Vector of R colors with length of 24.
#*		
getDefaultColors <- function()
{
    plot_colors <- c(552, 574, 645, 498, 450,  81,  26, 584, 
                     524, 472,  32,  57, 615, 635, 547, 254, 
                     100,  72, 630, 589,   8,  95, 568, 52);
	
    return (colors()[plot_colors]);
}


#* =======================================================================
#*		
#*    7)  fpFramePlot()
#*
#*    Plot frame distributions of ribosomal foot prints (bar plot)
#*
#*    Argument:
#*
#*        frame_file:  Character vector, name of the file with frame 
#*                     counts generated with fpFraming().
#*        by_column:   Logic, if the frame data arraned by column (row
#*                     is for each read length). If TRUE, fp_column
#*                     must be defined.
#*        fp_column:   Positive integer, column number for frame counts.
#*        read_len:    Positive integer vector, length of reads (alignments).
#*        title_text:  Character vector, text for title
#*
#*    Return:  None. Plot figure only
#*
fpFramePlot <- function(frame_file, by_column=TRUE, fp_column=c(2:4),
    read_len=c(25:30), legend_pos="topleft",
    title_text="Ribosomal fp Frame Distribution") {

    fp_frames <- read.delim(frame_file, header=TRUE);
    if(by_column == TRUE) {
        fp_frames <- t(fp_frames[, fp_column])
    }
    fp_frames <- fp_frames[, colnames(fp_frames) %in% read_len];
	
    frame_col <- c("red", "green",  "blue");
    barplot(fp_frames, beside=TRUE, col=frame_col,
        xlab="Read Length", main=title_text, cex.main=0.9);
    legend(legend_pos, legend=c("Frame 0", "Frame 1", "Frame 2"),
        fill=frame_col);
}


#*  ===========================================================================
#*
#*    8.  plotMetageneFrames()
#*
#*    Bar plot to metagene frame distribution.
#*
#*    Argument:
#*
#*        metagene_atStart: Character vector, name of the file which contains
#*                          frame name and frame counts for each metagene 
#*                          position relative to cds start position.	
#*        metagene_atStop:  Character vector, name of the file which contains
#*                          frame name and frame counts for each metagene 
#*                          position relative to cds stop position.
#*        min_5p:           Positive integer, minimum length of 5'end to start 
#*                          codon of cds.
#*        max_5p:           Positive integer, miximum length of 5'end to start 
#*                          codon of cds.
#*        min_3p:           Positive integer, minimum length of 3'end from cds 
#*                          start codon.
#*        max_3p:           Positive integer, miximum length of 3'end from cds 
#*                          start codon.
#*        frame_colors:     Character vector or R colors vector of length 3,
#*                          plot colors for reading frames.
#*        beside:           Logic, If FALSE, the columns of height are portrayed
#*                          as stacked bars, and if TRUE the columns will be 
#*                          portrayed as juxtaposed bars.
#*
#*    Return:  None. Plot figure only
#*	
plotMetageneFrames <- function(metagene_atStart=NULL, metagene_atStop=NULL,
    min_5p=-20, max_5p=200, min_3p=-200, max_3p=20,
    frame_colors=c("red", "green", "blue"), beside=TRUE) {

    if(is.null(metagene_atStart) | is.null(metagene_atStop))
        stop("Two input files are required.")
		
    to_start <- getRegionFrames(metagene_atStart, min_5p, max_5p);	
    to_stop  <- getRegionFrames(metagene_atStop,  min_3p, max_3p);			

    num_nas <- max(1, min(20, round(ncol(to_start) + ncol(to_stop) * 0.1)));
    NAs <- matrix(NA, ncol=num_nas, nrow=3);
    plot_data <- cbind(to_start, NAs, to_stop);
	
    bar_plot <- barplot(plot_data, beside=beside, col=frame_colors,
        ylim=c(0, max(plot_data, na.rm=TRUE)), 
        axisnames=FALSE, border=NA,
        space=c(0, 0.25), ylab="Number of Reads", 
        xlab="Base position relative to CDS");

    pretty5p <- getPrettyLabels(min_5p, max_5p);
    prettystart <- pretty5p;
    prettystart[prettystart == 0] <- "start";
    axis(side=1, at=as.vector(bar_plot)[(pretty5p - min_5p + 1)], 
        labels=prettystart, lwd=0, lwd.ticks=1);
                
    pretty3p <- getPrettyLabels(min_3p, max_3p);
    prettystop <- pretty3p;
    prettystop[prettystop == 0] <- "stop";
    axis(side=1, labels=prettystop, lwd=0, lwd.ticks=1,
        at = as.vector(bar_plot[, max(which(is.na(plot_data[1,])) + 
            1):ncol(bar_plot)])[pretty3p - min_3p + 1]);
}


#*  ===========================================================================
#*
#*    9.  getRegionFrames()
#*
#*    Extract a subset of metagene frames for bar plot.
#*
#*    Arguments:
#*
#*        frame_files:  Character vector, name of the file with metagene 
#*                      frame counts. The file is in tab-delimited foramt 
#*                      with 3 columns for counts, frame, and positions.
#*        from:         Integer, the first position for plot of metagene 
#*                      frames relative to cds start position.
#*        to:           Integer, the last position for plot of metagene 
#*                      frames relative to cds stop position.
#*
#*    Return:	
#*
#*        A data frame, subset of the metagene frame information and columns
#*        for position and rows for each reading frame
#*

getRegionFrames <- function(frame_file, from, to) {

    if(from > to) stop("'from' must be smaller than 'to'.");
    meta_frames <- read.table(frame_file, header=TRUE, 
						sep="\t", quote="");
    if(from >= max(meta_frames[,3], na.rm=TRUE))
        stop("Incorrect last position.");
		
    #*	Adjust 'from' and 'to' positions in order to start from 
    #*	the first frame0 read and stop at the last frame2 read
    #*	---------------------------------------------------------
    if(from %% 3 != 0) from <- from - from %% 3;
    if(to %% 3 != 2) to <- to + to %% 3;
	
    meta_frames <- meta_frames[meta_frames[,3] >= from,];
    meta_frames <- meta_frames[meta_frames[,3] <= to, ];
	
    #*	convert the table to three rows matrix. The first  
    #*	position must be the one for first frame 0. Note:
    #*	There is no guarantee that values in column 3 of 
    #*	the meta_frame (metagene position) are continuous
    #*	-------------------------------------------------
    positions <- from:to;
    frame_counts <- matrix(0, nrow=3, ncol=length(positions));

    frame0 <- meta_frames[which(meta_frames[,2] == 0), ];
    columns <- match(frame0[,3], positions);
    frame_counts[1,columns] <- frame0[,1];
	
    frame1 <- meta_frames[which(meta_frames[,2] == 1), ];
    columns <- match(frame1[,3], positions);
    frame_counts[2,columns] <- frame1[,1];
	
    frame2 <- meta_frames[which(meta_frames[,2] == 2), ];
    columns <- match(frame2[,3], positions);
    frame_counts[3,columns] <- frame2[,1];
	
    columns <- seq(from, to, 3);	
    frame_table <- matrix(0, nrow=3, ncol=length(columns));
    rownames(frame_table) <- c("frame0", "frame1", "frame2");
    colnames(frame_table) <- columns;	

    columns <- seq(1, ncol(frame_counts), 3);
    frame_table[1,] <- frame_counts[1, columns];
    frame_table[2,] <- frame_counts[2, columns+1];
    frame_table[3,] <- frame_counts[3, columns+2];
	
    return (frame_table);
}


#* ===========================================================================
#*
#*    10. getPrettyLabels()
#*
#*    Generate a pretty x axis label for metagene frame plot.
#*
#*    Arguments:
#*
#*        min_val:  Integer, minimum values of x axis range.
#*        max_val:  Integer, maximum values of x axis range.
#*
#*    Return: Integer vector, x axis label for metage frame plot.
#*

getPrettyLabels <- function(min_val, max_val) {

    pretty_lable <- pretty(min_val:max_val)
    pretty_lable <- pretty_lable[pretty_lable >= min_val];
    pretty_lable <- pretty_lable[pretty_lable <= max_val];
    pretty_lable <- unique(pretty_lable);
	
    return (pretty_lable);
}


#*  ===========================================================================
#*
#*    11.  metageneCountsPlot()
#*
#*    Arguments:
#*
#*        count_file:     Character vector, name (and path) of frame count file.
#*        from_position:  Integer, lestmost position of metagene to plot.
#*        to_position:    Integer, rightmost position of metagene to plot.
#*        codon:          Character vector, either "Start" or "Stop"
#*        x_interval:     Integer, length between tickmarks for x-axis
#*
#*    Return:  None. Plot figure only
#*

plotMetageneCounts <- function(count_file, from_position=-50, 
    to_position=50, codon="Start", x_interval=10) {

    if(from_position > to_position)
        stop("From position cannot be greater than to position.");
		
    meta_counts <- read.table(count_file, header=TRUE, sep="\t", quote="");
    from_row <- which(meta_counts[,3] == from_position)		
    last_row <- which(meta_counts[,3] == to_position)	

    plot(meta_counts$counts[from_row:last_row], type="l",
        ylab="Number of Reads", xaxt="n",
        xlab=paste("Distance to", codon, "Codon"));
	
    from_label <- from_position/x_interval;
    to_label <- to_position/x_interval;
    x_range <- length(from_position:to_position);
	
    x_labels <- c(c(from_label:-1)*x_interval, codon, 
        c(1:to_label)*x_interval);
    axis(1, at=seq(0, x_range, x_interval), labels=x_labels);
}


#*  =======================================================================
#*
#*    12)  plotHeatmap()
#*
#*    Make heatmap plot with data matrix. Since the image size may be
#*    big, the plot will be saved to output file.
#*
#*
#*    Arguments:
#*
#*        plot_value:   Numeric matrix, values for heatmap plot.
#*        sample_name:  Character vector, column (sample) labels.
#*        gene_name:    Character vector, row (gene) labels.
#*        image_type:   Character vector, output image format, 
#*                      supported are "pdf", "tiff", and "png".
#*        image_width:  Positive integer, width of output image
#*                      in inches.
#*        is.log2:      Logic, is the data log2 transformed.
#*        scale_by:     Character vector, how the data is scaled,
#*                      either "row" or "column"	
#*
#*    Return:  None. Plot figure only
#*
plotHeatmap <- function(plot_value, sample_name, gene_name, 
    image_type="pdf", image_width=8, is.log2=FALSE, scale_by="row") {

    image_height <- nrow(plot_value)/8 + 3;
    if(is.log2 == FALSE) plot_value <- log2(plot_value + 1);
	
    image_type <- tolower(image_type);
    if(image_type == "pdf") {
        pdf_file <- paste0(sample_name, ".cluster.pdf");
        pdf(pdf_file, width=image_width, height=image_height);
    } else if(image_type == "tiff") {
        tiff_file <- paste0(sample_name, ".cluster.tiff");
        tiff(tiff_file, width=image_width, height=image_height);
    } else if ((image_type == "png")) {
        png_file <- paste0(sample_name, ".cluster.png");
        png(png_file, width=image_width, height=image_height);
    } else { stop("The image type is not supported."); }
					
    heatmap.2(plot_value, scale=scale_by, trace="none", 
        col=bluered, density.info="none", srtCol=90,
        margins=c(12, 12), cexCol=1, cexRow=1,
        labRow=gene_name,
        key.title=NA, key.xlab=NA, key.ylab=NA,
        lhei=c(2, image_height-3, 0.75),  
        lwid=c(4,8,0.25),
        lmat=rbind(c(0, 3, 0), c(2, 1, 0), c(0, 4, 0))
    );
	
    dev.off();
}


#*  =======================================================================
#*
#*    13)  plotCorrelationHeatmap()
#*
#*    Plot sample correlation image with blue and red colors. Output
#*    will be save to file.
#*
#*    Arguments:
#*
#*        dds:         A DESeqDataSet object with size factor calculated 
#*                     or basepairs column presented in mcols.
#*        normalize:   Character vector, method of normalization either 
#*                     "size", "fpm", "fpkm", or NULL.
#*        image_name:  Character vector, output image file name.
#*        image_type:  Character vector, output image format, either
#*                     "pdf", "tiff", or "png".
#*        image_size:  Positive integer, image height.
#*
#*    Return:  None. Plot figure only
#*
plotCorrelationHeatmap <- function(dds, 
    normalize=c("size", "fpm", "fpkm"),
    image_name, image_type="pdf", image_size=12) {

    count_data <- getCountsData(dds, normalize);
    count_data <- log2(count_data + 1);

    sample_corr <- cor(count_data, 
        use=if (any(is.na(count_data))) "pairwise.complete.obs" else "all");
        rownames(sample_corr) <-colnames(count_data);
        colnames(sample_corr) <-colnames(count_data);
		
    plotRedBlueCorrelationImage(sample_corr, image_name, 
        image_type, image_size);			
}


#*  =======================================================================
#*
#*    14)  plotRedBlueCorrelationImage()
#*
#*    Plot sample correlation image with blue and red colors
#*
#*    Arguments:
#*
#*        corr_data:    Numeric matrix, pairwise correlation coefficients
#*                      of samples.
#*        image_name:   Character vector, output image file name.
#*        image_type:   Character vector, output image format, either
#*                      "pdf", "tiff", or "png".
#*        image_width:  Positive integer, image width.
#*
#*    Return:  None. Plot figure only
#*
plotRedBlueCorrelationImage <- function(corr_data, image_name, 
    image_type="pdf", image_width=12) {
		
    #  RGB colors for heatmap and color scale
    #  =====================================================
	
    RedRamp <- rgb(seq(1, 1, length=256),  		
                    seq(0, 1, length=256),  		
                    seq(0, 1, length=256)) ; 		

    BlueRamp <- rgb(seq(0, 1, length=256),  		
                    seq(0, 1, length=256),  		
                    seq(1, 1, length=256));  		
		
    ColorRamp   <- cbind(BlueRamp, rev(RedRamp));
    ColorLevels <- seq(min(corr_data), max(corr_data), 
                        length=length(ColorRamp));

    #*  PDF file width and height. Plot 5 samples
    #*  per inch for width and height
    #*  ============================================
	
    col_label <- colnames(corr_data);
    row_label <- rownames(corr_data);

    image_type <- tolower(image_type);
    if(image_type == "pdf") {
        pdf_file <- paste0(image_name, ".correlation_map.pdf");
        pdf(pdf_file, width=image_width, height=image_width);
    } else if(image_type == "tiff") {
        tiff_file <- paste0(image_name, ".correlation_map.tiff");
        tiff(tiff_file, width=image_width, height=image_width);
    } else if ((image_type == "png")) {
        png_file <- paste0(image_name, ".correlation_map.png");
        png(png_file, width=image_width, height=image_width);
    } else { stop("The image type is not supported."); }
			
    layout(matrix(data=c(1,2), nrow=1, ncol=2), 
        widths=c(9,1), heights=image_width);

    par(mai=c(2,2,2,2));
    image(1:length(corr_data[1,]),  1:length(corr_data[,1]),  
        zlim=c(min(corr_data), max(corr_data)),
        corr_data, col=ColorRamp,  
        xlab="", ylab="",  axes=FALSE,
        main=paste(image_name, "Sample Correlation"));
	
    for(label in 1:ncol(corr_data)) {
        axis(1, label, col_label[label], las=2, cex=0.5);
        axis(2, label, row_label[label], las=2, cex=0.5);
    }

    par(mai=c(4,0.1,4,0.75));
    image(2, ColorLevels, matrix(data=ColorLevels, 
        ncol=length(ColorLevels), nrow=1),
        col=ColorRamp, xlab="", ylab="", xaxt="n");

    dev.off();
}


#*  =======================================================================
#*
#*    15)  plotCorrelationMatrix()
#*
#*    This functin plot correlation matrix with corrplot() provided by
#*    corrplot package. For advanced R users, it is suggested to use 
#*    the corrplot package directly.
#*		
#*    The plot data could have values from one or different DESeqDataset 
#*    objects. For example, the normalized counts of Riboseq fp counts
#*    and RNAseq fp counts from one DESeqDataSet object, or normalized 
#*    fp counts of uORF regions and cds regions, cds length, 5'UTR length,
#*    TE change, ribosomal fp changes, mRNA fp changes. This should be
#*    generated and extracted from outside of this function.
#*			
#*    Arguments:
#*		
#*        plot_data:    A numeric matrix, NA is allowed but correlation 
#*                      matrix will contains NAs if NA exist.  
#*        cor_method:   Character vector, function name for correlation 
#*                      calculation, either "rcorr" (default) or "cor".
#*        cor_type:     Character vector, correlation type, either "spearman"
#*                      or "pearson". It is equal to the "type" argument for 
#*                      rcorr() or the "method" argument for cor().
#*        shape_type:   Character vector, shape to represent the correlation 
#*                      coefficient on the plot image, it is equal to the 
#*                       "method" argument for corrplot(), one of "circle",  
#*                      "square", "ellipse", "number", "shade",  "color", or 
#*                      "pie".
#*        p_threshold:  Numeric, the threshold for significant level
#*
#*    Return:  None. Plot figure only
#*
plotCorrelationMatrix <- function(plot_data, 
    cor_method="rcorr", cor_type="spearman", 
    shape_type="ellipse", p_threshold=0.01) {

    if(is.matrix(plot_data) == FALSE)
        stop("Plot data must be in a matrix.")
	
    if(cor_method == "cor") {
        RRR <- cor(plot_data, method=cor_type);
        corrplot(RRR, main="", method=shape_type);
    } else {
        cor_5 <- rcorr(as.matrix(plot_data), type="spearman")
        M <- cor_5$r;
        p_mat <- cor_5$P;
        corrplot(M, type="upper", order="hclust", 
            method=shape_type, p.mat=p_mat, 
            sig.level=p_threshold);
    }
}


#* =======================================================================
#*
#*    16)  plotMultiMetageneCounts()
#*		
#*    Plot reads counts along metagene positions for multiple samples
#*
#*    Arguments:
#*
#*        count_files:    Character vector, names (and path) of
#*                        metagene reads count files.
#*        from_position:  Integer, start position of metagene to plot. 
#*        to_position:    Integer, stop position of metagene to plot 
#*        codon:          Character vector, codon name, either "Start"
#*                        or "Stop" for x-axis labeling. 
#*        x_interval:     Positive integer, interval for x-axis lables
#*        line_colors:    Character vector of R color names for lines.
#*
#*    Return:  None. Plot figure only
#*

plotMultiMetageneCounts <- function(count_files, from_position=-50, 
    to_position=50, codon="Start", x_interval=10,
    line_colors=c("red", "blue", "green")) {

    if(from_position > to_position)
        stop("From position cannot be greater than to position.");
		
    if(length(count_files) > length(line_colors)) 
        line_colors <- rainbow(length(count_files));
	
    if(length(count_files) > 1 )
    count_matrix <- getMetageCountsMatrix(count_files, 
        from_position, to_position);
	
    for(a_col in seq_along(count_files)) {
        if(a_col == 1) {
            plot(count_matrix[, a_col], type="l", 
                xaxt="n", col=line_colors[a_col],
                xlab=paste("Distance to", codon, "Codon"),
                ylab="Number of Reads" );
				
            from_label <- from_position/x_interval;
            to_label <- to_position/x_interval;
            x_range <- length(from_position:to_position);

            x_labels <- c(c(from_label:-1)*x_interval, codon, 
                c(1:to_label)*x_interval);
            axis(1, at=seq(0, x_range, x_interval), labels=x_labels);			
				
        } else {
            lines(x=1:x_range, y=count_matrix[, a_col], 
                col=line_colors[a_col]);
        }
    }	
}


#*  =======================================================================
#*
#*    17)  getMetageCountsMatrix()
#*
#*    Construct count matrix from metagene frame tables 
#*
#*    Arguments:
#*
#*        count_files:    Character vector, names (and path) of
#*                        metagene reads count files.
#*        from_position:  Integer, start position of metagene to plot. 
#*        to_position:    Integer, stop position of metagene to plot 
#*
#*    Return:
#*
#*        Data matrix with read counts along metagene position from
#*        multiple samples. Columns in matrix is sorted by
#*        maximum counts of each sample in decreasing order.
#*

getMetageCountsMatrix <- function(count_files, from_position, 
    to_position) {

    if(length(count_files) == 1) 
        stop("Count files must be more than one.");

    meta_counts <- matrix(0, ncol=length(count_files), 
        nrow=length(from_position:to_position));
    max_counts <- rep(0, length(count_files));
    for(a_file in seq_along(count_files)) {
        in_counts <- read.table(count_files[a_file], header=TRUE, 
            sep="\t",quote="");
		
        from_row <- which(in_counts[,3] == from_position);
        last_row <- which(in_counts[,3] == to_position);	
		
        if(length(from_row) == 0 || length(last_row) == 0)
            stop("Input file has different metagene locations.");
		
        in_counts <- in_counts$counts[from_row:last_row]
        meta_counts[,a_file] <- in_counts;
        max_counts[a_file] <- max(in_counts);
    }
	
    col_order <- order(max_counts, decreasing=TRUE);
    meta_counts <- meta_counts[, col_order]
}

#*
#*  End of RiboProfiling_Plot_Util_Source.R
#*  __________________________________________________________________
#*  <Plot_Util><Plot_Util><Plot_Util><Plot_Util><Plot_Util><Plot_Util>

	
