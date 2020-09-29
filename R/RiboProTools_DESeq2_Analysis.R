#*
#*		Functions for DESeq2 analysis with ribosomal profiling data
#*
#*		Function implemented:
#*
#*		1.		getCountMatrixFromFiles()
#*		2.		getCountMatrixFromTable()
#*		3.		filterCountMatrix()
#*		4.		getAnnotationInfo()
#*		5.		getDESeqDataSet()
#*		6.		extractEfficiencyChange()
#*		7.		extractTranscriptionChange()
#*		8.		quickDESeq2Test()
#*		9.		getUORFCountMatrixFromTwoFileSets()
#*		10.		getUorfCountMatrixInFujunFormat()
#*		11.		runDESeq()
#*		12.		getUORFCountMatrixFromOneFileSet()
#*
#*		Last revised on April 18, 2018
#*		____________________________________________________________________
#*		<DESeq2 Analysis><DESeq2 Analysis><DESeq2 Analysis><DESeq2 Analysis>




#*	=======================================================================
#*
#*		1.	getCountTableFromRawCountFiles()
#*		
#*		Read each raw count file and merge them as one matrix
#*
#*		Arguments:
#*
#*			directory_name:	Character vector, the directory name  
#*								where raw count files can be found.
#*			file_name_pattern:	Character vector, the common pattern of 
#*								raw count files.
#*			count_column:		Non-negative integer, the column of raw 
#*								count value.
#*			rowname_column:	Non-negative integer, the column in raw  
#*								count file for row names. 0 for raw count 
#*								files saved by write.table() with 
#*								row.names=TRUE.
#*			has.header:			Logical, if the columns has headers.
#*		
#*		Return:
#*
#*			A data matrix with row names
#*
#*
getCountMatrixFromFiles <- function(directory_name, file_name_pattern, 
				count_column, rowname_column=0, has.header=TRUE)
{
	count_files <- list.files(pattern=file_name_pattern,
						path=directory_name, full.names=TRUE);

	count_table <- NULL;
	for(a_file in 1:length(count_files)) 
	{
		count_in <- read.table(count_files[a_file], header=has.header, 
						sep="\t", quote="");
		count_table <- cbind(count_table, count_in[,count_column]);
	}
	
	colnames(count_table) <- basename(count_files);
	
	if(rowname_column == 0) {
		rownames(count_table) <- rownames(count_in);
	} else {
		rownames(count_table) <- count_in[, rowname_column];
	}
	
	return (count_table);
}


#'	=======================================================================
#*
#*		2.	getCountMatrixFromTable()
#* 	
#*		In general, Ribosomal profiling data processing will generate 
#*		a count table that contains both Riboseq and RNAseq counts for 
#*		same contol and treatment samples. This function will extract
#*		a matrix from count table by columns if there are more than two 
#*		treatment groups , or convert the count table to a matrix, and
#*		finally reorder the columns to put control samples first then 
#*		treatment samples.
#*
#*		Arguments:
#*
#*			count_table:	A data frame or matrix with raw counts from 
#*							both Riboseq and RNAseq data for same samples.
#*							Row names must be gene names or gene IDs.
#*			ribo_control: 	Positive integer, columns for ribosomal control
#*							samples.
#*			ribo_treatmen: 	Positive integer, columns for ribosomal 
#*							treatment samples.
#*			mRNA_control: 	Positive integer, columns for mRNA control 
#*							samples.
#*			mRNA_treatment:	Positive integer, columns for mRNA treatment 
#*							samples.
#*		
#*		Return:		
#*
#*			A matrix with columns in the order of:
#*
#*				mRNA control samples,
#*				mRNA treatment samples, 
#*				ribosomal control samples, 
#*				ribosomal treatment samples
#*
#*
getCountMatrixFromTable <- function(count_table=NULL, 
			ribo_control=NULL, ribo_treatment=NULL,
			mRNA_control=NULL, mRNA_treatment=NULL)
{
	#*	All arguments must be provided
	#*	------------------------------
	if(is.null(count_table) || is.null(ribo_control) || 
		is.null(ribo_treatment) || is.null(mRNA_control) || 
		is.null(mRNA_treatment))
		stop("Missing one or more function argument.");
	
	#*	count_table must be data frame or matrix with
	#*	numeric values
	#*	-----------------------------------------------------
	if(!is.matrix(count_table) && !is.data.frame(count_table) )
		stop("count_table must be a data frame or data matrix.");
	if(!is.numeric(as.matrix(count_table)))
		stop("count_table must be numeric.")
			
	all_sample_columns <- c(ribo_control, ribo_treatment,
			mRNA_control, mRNA_treatment);

	#*	column numbers for each group must be unique
	#*	--------------------------------------------	
	if(length(unique(all_sample_columns)) < length(all_sample_columns))
		stop("One or more columns is defined in two groups.");
		
	#*	column numbers cannot be 0 or negative
	#*	--------------------------------------
	if(length(which(all_sample_columns < 1)) > 0 )
		stop("Column number cannot be 0 or negative.");
		
	#*	column numbers must be in range of count_table
	#*	----------------------------------------------
	if(max(all_sample_columns) >  ncol(count_table))
		stop("Sample columns are out of range.");
	
	count_matrix <- cbind(count_table[,mRNA_control], 
						  count_table[,mRNA_treatment], 
						  count_table[,ribo_control], 
						  count_table[,ribo_treatment]);
	rownames(count_matrix) <- rownames(count_table);
	
	return (as.matrix(count_matrix));
}


#'	=======================================================================
#*
#*		3.	filterCountMatrix()
#*		
#*		Filter the count matrix by mRNA mean or ribo mean to remove 
#*		the rows with low count of mRNA reads or ribo footprints
#*	
#*		Arguments:
#*
#*			count_matrix: A data matrix
#*
#*			mRNA_col:	Positive integer vector, the columns for
#*						mRNA reads counts in count_matrix.
#*			ribo_col:	Positive integer vector, the columns for 
#*			 			ribosomal footprints counts in count_matrix.
#*			mRNA_level:	Non negative numeric, mRNA values used to filter 
#*						out of data based on row means of mRNA counts.
#*			ribo_level:	Non negative numeric, ribo values used to filter 
#*						out of data based on row means of ribo counts.
#*
#*		Return:
#*
#*			Same data matrix as input but low count rows removed.
#*
#*
filterCountMatrix <- function(count_matrix=NULL, 
					mRNA_col=NULL, ribo_col=NULL, 
					mRNA_level=10, ribo_level=0)
{
	#*	Arguments without default values must be provided
	#*	--------------------------------------------------------------
	if(is.null(count_matrix) || is.null(mRNA_col) || is.null(ribo_col)) 
		stop("Missing count matrix, or mRNA_col or ribo_col.");
	
	#*	All arguments must be numeric
	#*	--------------------------------------------------
	if(!is.numeric(count_matrix) || 
		!is.numeric(mRNA_col) ||  !is.numeric(ribo_col) ||  
		!is.numeric(mRNA_level) || !is.numeric(ribo_level)) 
		stop("All srguments must be numeric.");
	
	#*	No negative number for mRNA_col or ribo_col, column 
	#*	numbers must be in range of count_matrix dimenssion,
	#*	and mRNA_col or ribo_col cannot have overlap
	#*	-------------------------------------------------------
	col_numbers <- c(mRNA_col, ribo_col);
	if(min(col_numbers) <= 0)
		stop("Column numbers cannot be 0 or negative.");
	if(max(col_numbers) > ncol(count_matrix))
		stop("Column number is out of count_matrix dimension.");
	if(length(col_numbers) != length(unique(col_numbers)))
		stop("mRNA_col and ribo_col have overlap.")
	
	#*	mRNA_level and ribo_level must be non-negative
	#*	----------------------------------------------
	if(mRNA_level < 0 || ribo_level < 0)
		stop("Thresholds must be 0 or greater.")
	if(mRNA_level > 0) {
		mRNA_mean <- rowMeans(count_matrix[, mRNA_col]);
		count_matrix <- count_matrix[mRNA_mean >= mRNA_level,];
	}
	if(ribo_level > 0) {
		ribo_mean <- rowMeans(count_matrix[, ribo_col]);
		count_matrix <- count_matrix[ribo_mean >= ribo_level,];
	}

	return (count_matrix);
}


#*	=======================================================================
#*
#*		4.	getAnnotationInfo()
#*
#*		Extract annotation data from a database file for selected genes .
#*
#*		Argument:
#*
#*			annotation:	A data frame of 3 columns for gene ID, gene  
#*							name, and gene description. 
#*
#*			gene_list:		Character vector, list of gene ID or gene 
#*							names for which annotation info will be 
#*							extracted.
#*
#*		Return:
#*
#*			A data frame of columns 3 for gene ID, name and description
#*
#*
getAnnotationInfo <- function(annotation, gene_list, id_column=1,   
		name_column=2, description=3)
{
	if(ncol(annotation) < 3)
		stop("Annotation does not have enough columns.");
	
	gene_id <- as.character(annotation[, id_column]);
	gene_rows <- which(gene_list %in% gene_id);
	if(length(gene_rows) != length(gene_list))
		stop("Cannot match all genes.")
		
	rows <- match(gene_list, as.character(annotation[, id_column]));
	columns <- c(id_column, name_column, description);
	gene_info <- annotation[rows, columns];
	colnames(gene_info) <- c("gene_ID","gene_name", "description");

	return (gene_info);
}


#*	=======================================================================
#*
#*		5.	getDESeqDataSet()
#*
#*		Generate a DESeqDataSet from count matrix with fixed order of 
#*		sample columns, fixed colData, and full design model:
#*			 ~ genotype + condition + genotype:condition
#*
#*		Arguments:
#*
#*			count_matrix:		A matrix with columns in the order of:
#*								mRNA control samples, and mRNA treatment  
#*								samples, ribosomal control samples, and  
#*								ribosomal treatment samples.
#*			num_Ribo_wildtype:	Positive integer, total number of ribo 
#*								wildtype samples.
#*			num_Ribo_mutant: 	Positive integer, total number of ribo 
#*								mutant samples.
#*			num_mRNA_wildtype:	Positive integer, total number of mRNA 
#*								wildtype samples.
#*			num_mRNA_mutant: 	Positive integer, total number of mRNA 
#*								mutant samples.
#*			annotation_info:	A data frame with columns of 2 for gene 
#*								name and gene description, and with gene 
#*								ID as row names.
#*
#*		Return:
#*
#*			A DESeqDataSet object with interaction term design from two
#*			conditions and two groups, as well as gene descriptions.
#*	
#*
getDESeqDataSet <- function(count_matrix, 
			num_Ribo_wildtype, num_Ribo_mutant, 
			num_mRNA_wildtype, num_mRNA_mutant, 
			annotation_info=NULL)
{
	assay <- c(rep("mRNA", time=num_mRNA_wildtype + num_mRNA_mutant),
			  rep("Ribo", time=num_Ribo_wildtype + num_Ribo_mutant));
	genotype <- c(rep("wildtype", time=num_mRNA_wildtype), 
			  rep("mutant", time=num_mRNA_mutant), 
			  rep("wildtype", time=num_Ribo_wildtype), 
			  rep("mutant", time=num_Ribo_mutant) );

	experiment_info <- DataFrame(
		condition=factor(assay, levels=c("mRNA", "Ribo")),
		genotype=factor(genotype, levels=c("wildtype", "mutant")),
		row.names=colnames(count_matrix)
	);
	
	DESeqDataSet <- DESeq2::DESeqDataSetFromMatrix(
				countData=count_matrix, colData=experiment_info, 
				design= ~ genotype + condition + genotype:condition);
	if(is.null(annotation_info) == FALSE)
		mcols(DESeqDataSet) <- annotation_info;
	
	return (DESeqDataSet);
}


#*	=======================================================================
#*
#*		6.	extractEfficiencychange()
#*
#*		Extract translational efficiency (TE) and TE changes with
#*		interaction term design
#*
#*		Arguments:
#*
#*			dds_object:			A DESeq object returned from DESeq(dds) or
#*								nbinomWaldTest(dds). The design model must be
#*								~ genotype + condition + genotype:condition.
#*			control_name:		Character vectors, name of control group.
#*			mutant_name:		Character vectors, name of mutant group.
#*			efficiency_type:	Character vectors, name of efficiency
#*								type, such as "TE" or "RRO"
#*			meta_cols:			Positive integer, columns of meta-data to be 
#*								attached to outputs, set 0 for nothing to 
#*								attach, NULL to output all columns.
#*		Return:
#*
#*			A DESeqDataSet with new elements added. 
#*
#*		Note: 	All results are save in file and returned object
#*			 	is for saving and debugging purposes.
#*
#*/
#'
extractEfficiencyChange <- function(dds_object, control_name, 
			mutant_name, efficiency_type="TE", meta_cols=c(1:3))
{
	TE_pattern  <- paste0("_", efficiency_type, "_with_");
	TEC_pattern <- paste0(efficiency_type, "_change_");
	dds <- dds_object;
	
	if(is.null(meta_cols)) meta_cols <- 1:ncol(mcols(dds));
	
	#*		The translational efficiency of control sample(s)
	#*		(condition effect for genotype I), use 
	#*
	#*		results(dds, contrast=c("condition","Ribo","mRNA"))  # or 
	#*		results(dds, name="condition_Ribo_vs_mRNA")

	control_TE <- DESeq2::results(dds, contrast=c("condition","Ribo","mRNA"));
	result2save_control <- data.frame(control_TE, mcols(dds)[, meta_cols]);

	ctl_out <- paste0(control_name, TE_pattern, mutant_name, ".txt");
	write.table(result2save_control, file=ctl_out, sep="\t", 
				quote=FALSE, row.names=TRUE, col.names=TRUE);


	#* 	The translational efficiency of mutant sample(s) 
	#*	(condition effect for genotype II)
	#*	Per example of results() documentation:
	#*	results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))

	mutant_TE <- DESeq2::results(dds, 
		list(c(resultsNames(dds)[3], resultsNames(dds)[4])));
	result2save_mutant <- data.frame(mutant_TE, mcols(dds)[,meta_cols]);

	mutant_out <- paste0(mutant_name, TE_pattern, control_name, ".txt");
	write.table(result2save_mutant, file=mutant_out, sep="\t", 
		quote=FALSE, row.names=TRUE, col.names=TRUE);

	#*	The TEC (translational efficiency change) of mutant sample(s):
	#* 	the interaction term, answering: is the condition effect *different* 
	#* 	across genotypes? 
	#*	--------------------------------------------------------------------
	#* 		results(dds, name="genotypeII.conditionB") 

	TEC_mutant <- DESeq2::results(dds, name=resultsNames(dds)[4]);
	result2save_TEC <- data.frame(TEC_mutant, mcols(dds)[,meta_cols]);
	
	TEC_out <- paste0(TEC_pattern, mutant_name, "_vs_", control_name, ".txt");
	write.table(result2save_TEC, file=TEC_out, sep="\t",	
		 quote=FALSE, row.names=TRUE, col.names=TRUE);

	return (dds);
}


#*	=======================================================================
#*
#*		7.	extractTranscriptionChange()
#*
#*		Extract mRNA and Ribo changes with group variables.
#*
#*		Arguments:
#*
#*			dds_object:		A DESeq object returned from running of 
#*							DESeq(dds) or nbinomWaldTest(dds). The 
#*		 				 	design model must be ~ group.	
#*			control_name:	Character vectors, name of control group.
#*			mutant_name:	Character vectors, name of mutant group.
#*	
#*		Return:
#*
#*			DESeqDataSet with new elements added. 
#*
#*
#*		Note: 	All results are save in file and returned object
#*			 	is for saving and debugging purposes.
#*
#*/
extractTranscriptionChange <- function(dds_object, 
			control_name, mutant_name, change_one="mRNA", 
			change_two="Ribo", meta_cols=c(1:3))
{
	pattern_one <- paste0("_", change_one, "_change_vs_");
	pattern_two <- paste0("_", change_two, "_change_vs_");

	dds <- dds_object;
	pair_names <- unique(as.character(dds$group));
	
	if(is.null(meta_cols)) meta_cols <- 1:ncol(mcols(dds));
	
	mRNA_change <- DESeq2::results(dds, 
		contrast=c("group", pair_names[2], pair_names[1]));
	result2save_mRNA <- data.frame(mRNA_change, mcols(dds)[,meta_cols]);
	
	mRNA_out <- paste0(mutant_name, pattern_one, control_name,".txt");
	write.table(result2save_mRNA, file=mRNA_out, sep="\t",	quote=FALSE, 
					row.names=TRUE, col.names=TRUE);
		
	Ribo_change <- DESeq2::results(dds,
		contrast=c("group", pair_names[4], pair_names[3]));
	result2save_Ribo <- data.frame(Ribo_change, mcols(dds)[,meta_cols]);
	
	Ribo_out <- paste0(mutant_name, pattern_two, control_name,".txt");
	write.table(result2save_Ribo, file=Ribo_out, sep="\t",	quote=FALSE, 
			row.names=TRUE, col.names=TRUE);

	return (dds);
}


#*	=======================================================================
#*	
#*		8.	quickDESeq2Test()
#*
#*		Perform DESeq analysis with one call
#*
#*		Arguments:
#*
#*			count_table:	 	Numeric matrix with mRNA counts and Ribo 
#*								Fp counts.
#*			Ribo_wildtype:		Positive integer vectors, column numbers of 
#*							 	counts for Ribo_wildtype in count matrix.
#*			Ribo_mutants: 	 	Positive integer vectors, column numbers of
#*							 	counts for Ribo_mutants in count matrix.
#*			mRNA_wildtype:	 	Positive integer vectors, column numbers of 
#*							 	counts for mRNA_wildtype in count matrix.
#*			mRNA_mutants:	 	Positive integer vectors, column numbers of 
#*							 	counts for mRNA_mutants in count matrix.
#*			control_name:	 	Character vector, names of control group. 
#*			mutant_genotype: 	Character vector, names of mutant group.
#*			mRNA_level:		 	Numeric, thresholds to filter out the matrix
#*							 	based on row means of mRNA samples.
#*			Ribo_level:		 	Numeric, thresholds to filter out the matrix
#*							 	based on row means of Ribo samples.			
#*			annotation_file:	Character vector, name of gene annotation file
#*	
#*		
#*		Return:	None. All outputs will be saved in files.
#*		
#*
quickDESeq2Test <- function(count_table, 
			Ribo_wildtype, Ribo_mutants, 
			mRNA_wildtype, mRNA_mutants, 
			control_name, mutant_name,
			mRNA_level, Ribo_level,			
			annotation_file)
{
	count_matrix <- getCountMatrixFromTable(count_table, 
			Ribo_wildtype, Ribo_mutants, 
			mRNA_wildtype, mRNA_mutants);

	mRNA_col <- 1:sum(length(mRNA_wildtype) + length(mRNA_mutants));
	ribo_col <- mRNA_col + sum(length(Ribo_wildtype) + length(Ribo_mutants));
	
	count_matrix <- filterCountMatrix(count_matrix, 
			mRNA_col, ribo_col,
			mRNA_level, Ribo_level);

	gene_list <- rownames(count_matrix);
	gene_info <- getAnnotationInfo(annotation_file, gene_list);
	
	num_Ribo_wildtype <- length(Ribo_wildtype);
	num_Ribo_mutant   <- length(Ribo_mutants);
	num_mRNA_wildtype <- length(mRNA_wildtype);
	num_mRNA_mutant   <- length(mRNA_mutants);

	dds <- getDESeqDataSet(count_matrix, num_Ribo_wildtype, 
			num_Ribo_mutant, num_mRNA_wildtype, num_mRNA_mutant, 
			gene_info);

	dds_TE <- extractEfficiencyChange(dds, control_name, mutant_name);
	TE_file <- paste0(mutant_name, "_vs_", control_name, "_TE_change.RData");
	save(dds_TE, file=TE_file);
	
	dds_RC <- extractTranscriptionChange(dds, control_name, mutant_name);
	RC_file <- paste0(mutant_name, "_vs_", control_name, "_riboChange.RData")
	save(dds_RC, file=RC_file);
	
	normalized_counts <- DESeq2::counts(dds_TE, normalized=TRUE);
	write.table(normalized_counts, sep="\t", quote=FALSE,
		row.names=TRUE, col.names=TRUE,
		file=paste0(mutant_name,"_normalized_counts_from_DESeq2.txt"));
}


#*		=====================================================================
#*
#*		9.	getUORFCountMatrixFromTwoFileSets()
#*
#*		Generate count table from single count files for DESeq2 analysis 
#*		with default normalization or normalization by subsets. Count 
#*		files must include both uorf count files and gene cds count files.
#*		Unique file name patterns should be contained in uorf count files
#*		and gene cds count files for distinguishing them each other. Also
#*		row names of uorf counts must contain gene name and uorf name, and
#*		a common pattern between gene name and uorf name. Each gene name
#*		in row names of uorf counts must be in row names of cds counts.
#*					
#*		Argument:
#*
#*			file_directory:		Character vector, directory where single  
#*								count files are stored.
#*			ourf_file_pattern:	Character vectors, pattern of ourf count 
#*								file names.
#*			cds_file_pattern:	Character vectors, pattern of cds count 
#*								file names.
#*			count_column: 		Possitive interger greater than 0. The  
#*								column of raw count in count file.
#*			rowname_column:		Possitive interger, which column will be  
#*								used as row names. Set to 0 if count files 
#*								are output from R write.table() with 
#*								row.names=TRUE. 
#*			has.header:			Logic, if the input files has headers.
#*			seperator:			Character vector, unique pattern to separate 
#*								gene name from uorf names. It CANNOT be any 
#*								character in part of gene name or uorf id 
#*								in uorf names.
#*		Return:
#*
#*			A data matrix with raw footprints counts. The left half of the 
#*			matrix is fp counts of uorf regions for each sample and the 
#*			right half is fp counts of matched gene cds regions.				
#*
#*
getUORFCountMatrixFromTwoFileSets <- function(file_directory, 
			uorf_file_pattern="uorf_conserved_fpcounts.txt", 
			cds_file_pattern="cds_for_uorf_fpcounts.txt",
			count_column=2, rowname_column=0, 
			has.header=TRUE, seperator="\\.")
{
	#*	fp counts based on uorf regions
	#*	-----------------------------------------------------
	uorf_counts <- getCountMatrixFromFiles(file_directory, 
				uorf_file_pattern, count_column=2, 
				rowname_column=0, has.header=TRUE);
	colnames(uorf_counts) <- sub(uorf_file_pattern, "_uorf",
							colnames(uorf_counts));

	#*	fp counts based on transcript cds regions
	#*	-----------------------------------------------------
	cds_counts  <- getCountMatrixFromFiles(file_directory, 
							cds_file_pattern, count_column, 
							rowname_column, has.header);
	colnames(cds_counts) <- sub(cds_file_pattern, "_cds",
							colnames(cds_counts));
	
	#*	match rows by gene name and merge two matrix as one
	#*	-------------------------------------------------------
	cds_gene <- rownames(cds_counts);
	uorf_gene <- sub(paste0(seperator, ".{,1}"), "", 
						rownames(uorf_counts));

	cds_rows <- match(uorf_gene, cds_gene);
	if(length(is.na(cds_rows)) > 0)
		message("Some uorf do not match any gene.")
	count_matrix <- cbind(uorf_counts, cds_counts[cds_rows]);
	
	return (count_matrix);
}


#*		=====================================================================
#*
#*		10.	getUorfCountMatrixInFujunFormat()
#*
#*		Generate count table from single count files with Fujun*s method
#*		for DESeq2 analysis. Each count file contains both uorf fp counts
#*		and gene cds fp counts, usually based on a bed file containing
#*		both uorf and cds definition. Each row name for uorf rows must 
#*		have a gene name and uorf identifier, and the gene names in uorf
#*		rows must be in row names of cds counts.
#*
#*		The final count matrix will have four parts:
#*
#*		1)	Top left:		fp counts of uorf region for each sample
#*		2)	Top right:		fp counts of matched gene cds region
#*		3)	Bottom left:	fp counts of all gene cds region
#*		4)	Bottom right:	fp counts of all gene cds region
#*
#*		Argument:
#*
#*			file_directory:	Character vector, directory where single 
#*								count files are stored.					
#*			uorf_file_pattern:	Character vectors, pattern in names of uorf  
#*								fp count files.
#*			mRNA_file_pattern:	Character vectors, pattern in names of mRNA 
#*								count files. If null, relavent part will be 
#*								extracted from cds part of uorf counts.
#*			count_column: 		Possitive interger greater than 0. the column 
#*								of raw count in count file.
#*			rowname_column: 	Possitive interger, which column will be used 
#*								as row names. Set to 0 if count files are 
#*								output from R write.table() with row.names=TRUE.
#*			has.header: 		Logic, if both count files have column headers.
#*			separator:			Character vector, unique pattern in uorf names 
#*								used to extract gene names. uorf names must 
#*								follow the format of gene_name + separator + 
#*								unique uorf id. The seperator cannot be any
#*								pattern existing in any gene names.
#*
#*		Return:
#*		
#*			A data frame with contents described above.
#*
#*
getUORFCountMatrixInFujunFormat <- function(file_directory, 
			uorf_file_pattern="uorf_conserved_fpcounts.txt", 
			cds_file_pattern="cds_for_uorf_fpcounts.txt",
			count_column=2, rowname_column=0, has.header=TRUE)
{
	#*	fp counts based on uorf regions
	#*	-----------------------------------------------------
	uorf_counts <- getCountMatrixFromFiles(file_directory, 
				uorf_file_pattern, count_column=2, 
				rowname_column=0, has.header=TRUE);
	colnames(uorf_counts) <- sub(uorf_file_pattern, "_uorf",
							colnames(uorf_counts));

	#*	fp counts based on transcript cds regions
	#*	-----------------------------------------------------
	cds_counts  <- getCountMatrixFromFiles(file_directory, 
							cds_file_pattern, count_column, 
							rowname_column, has.header);
	colnames(cds_counts) <- sub(cds_file_pattern, "_cds",
							colnames(cds_counts));
	
	#*	merge two matrix as one
	#*	-------------------------------------------------------
	count_matrix <- cbind(uorf_counts, cds_counts);

	return (count_matrix);
}




#*	=======================================================================
#*	
#*		11.	runDESeq()
#*
#*		Get DESeq object from DESeq() with a DESeqDataSet returned 
#*		from getDESeqDataSet() or perform estimateDispersion() and  
#*		nbinomWaldTest() with the DESeqDataSet returned from the
#*		getDESeqDasetWithSubsetSizeFactors() function call.
#*			  
#*		Arguments:
#*	
#*			deseq_dataset:		A DESeq object returned from getDESeqDataSet() 
#*								or from getDESeqDasetWithSubsetSizeFactors().
#*			has.SizeFactors:	Logic, if the deseq_daaset has size factor.
#*			reset.design:		Logic, if need reset design model. The 
#*								default design model is with interactions:
#*								~ genotype + condition + genotype:condition 
#*								whihc is for translational efficiency test. 
#*								For tranascription test, it must be changed 
#*								to ~group
#*		Return:
#*
#*			A DESeq object with all test results
#*/

runDESeq <- function(deseq_dataset, has.SizeFactors=FALSE,
		reset.design=FALSE, fit_type="parametric")
{
	dds <- deseq_dataset;
	
	#*	new design by group variables
	#*	-----------------------------------------
	if(reset.design == TRUE) {
		group <- paste0(dds$genotype, dds$condition)
		dds$group <- factor(group);
		design(dds) <- ~ group;	
	}

	if(has.SizeFactors == FALSE) {
		dds <- DESeq2::DESeq(dds, fitType=fit_type);
	} else {
		dds <- DESeq2::estimateDispersions(dds, fitType=fit_type);
		dds <- DESeq2::nbinomWaldTest(dds);
	}
	
	return (dds);
}


#*	=======================================================================
#*	
#*		12.		getUORFCountMatrixFromOneFileSet()
#*
#*		Generate uorf count matrix from files which contains fp counts
#*		for uorf and cds based on one bed file. Genes of uorf must be
#*		included in both uorf rows and gene rows such as defined in 
#*		yeast_all_plus_uorf2721.bed file.
#*
#*		Arguments:
#*
#*			directory_name:		Character vector, the directory name where 
#*								raw count files can be found. 
#*			file_name_pattern:	Character vector, the common pattern of raw 
#*								count files.
#*			count_column:		Positive integer, column of raw count value.
#*			rowname_column:		Positive integer, the column in raw count 
#*								file for row names. 0 for raw count files 
#*								saved by write.table with row.names=TRUE.
#*			has.header:			Logical, if the raw count files have headers.
#*			gene_ID_patter:		Character vector, unique pattern of gene ID 
#*								used to distinguish gene ID from uorf ID.
#*			uorf_ID_patter:		Character vector, unique pattern of uorf ID 
#*								used to distinguish uorf ID from gene ID.
#*
#*
getUORCountMatrixFromOneFileSet <- function(directory_name, 
			uorf_file_pattern, mrna_file_pattern=NULL, 
			count_column=2, rowname_column=0, has.header=TRUE, 	
			uorf_ID_patter1="^.{4}-", uorf_ID_patter2=NULL)
{
	uorf_counts <- getCountMatrixFromFiles(directory_name, 
							uorf_file_pattern, count_column, 
							rowname_column, has.header);
	
	colnames(uorf_counts) <- paste0(colnames(uorf_counts), "_uorf");
	
	uorf_genes <- sub(uorf_ID_patter1, "", rownames(uorf_counts));
	if(!is.null(uorf_ID_patter2)) 
		uorf_genes <- sub(uorf_ID_patter2, "", uorf_genes);
	
	if(is.null(mrna_file_pattern)) {
		cds_rows <- match(uorf_genes, rownames(uorf_counts));
		cds_part <- uorf_counts[cds_rows,];
		colnames(cds_part) <- sub("_uorf$", "_cds", colnames(cds_part));
	} else {
		mrna_counts <- getCountMatrixFromFiles(directory_name, 
							mrna_file_pattern, count_column, 
							rowname_column, has.header);
		colnames(mrna_counts) <- paste0(colnames(mrna_counts), "_mrna");
		
		cds_rows <- match(uorf_genes, rownames(mrna_counts));
		cds_part <- mrna_counts[cds_rows,];
	}
	
	count_matrix <- cbind(uorf_counts, cds_part);

	return (count_matrix);
}


#*
#*	End of Ribosomal_Profiling_DESeq2_Analysis_Source.R
#*
#*	=========================================================================










