####
## Functions for reading in GTF files
## By: Olivia
## Last major edit: Mar 2 2021
####

# requires: the temp_unzip() function from source_temp_unzip.R


####
## Reading in GTF
read_gtf <- function(filename, data.table = F, keep.extra.columns = T){
  if (!file.exists(filename)) {
    stop("No such file: ", filename);
  }
  print(paste("Reading in", filename))
  if (endsWith(filename, "gz") | endsWith(filename, "zip") | endsWith(filename, "bzip2") | endsWith(filename, "xz")) {
    output <- suppressWarnings(temp_unzip(filename, fread, data.table = data.table))
  } else {
    output <- suppressWarnings(fread(filename, data.table = data.table))
  }
  print("GTF read. Now checking that it conforms to GTF structure...")
  # TEST FILE STRUCTURE:
  # should be at least 9 columns (if more, note how many extra there are for column naming later)
  num_extra_columns <- 0
  if (ncol(output) < 9) {
    stop("Fewer than 9 columns detected")
  }
  if (ncol(output) > 9) {
    num_extra_columns <- ncol(output) - 9
  }
  # all of column 1 should start with "chr"; 
  if (!all(startsWith(output[,1], "chr") == T)) {
    stop("Error in GTF: not all entries in column 1 (seqname) are chr--")
  }
  # column 3 should include c("gene", "transcript", "exon"); 
  if (!(any(output[,3] == "gene") & any(output[,3] == "transcript") & any(output[,3] == "exon"))) {
    stop("Error in GTF: missing either gene, transcript, or exon level information in column 3 (type)")
  }
  # columns 4 and 5 should be integers
  if (!(is.integer(output[,4]) & is.integer(output[,5]))) {
    stop("Error in GTF: columns 4 (start) and 5 (end) are not integers")
  }
  # column 4 should be <= column 5; 
  if (!all(output[,4] <= output[,5])) {
    stop("Error in GTF: column 4 (start) is not always <= column 5 (end)")
  }
  # column 7 should all be "+" or "-" ;
  if (!all(output[,7] %in% c("+", "-"))) {
    stop("Error in GTF: column 7 is not defining strand")
  }
  # column 9 contains a bunch of stuff, but a basic check is "gene_id" should always be there
  if (!all(grepl("gene_id", output[,9]))) {
    stop("Error in GTF: column 9 (attribute) is missing gene_id, the minimum info")
  }
  print("Tests passed!")
  # make informative column names
  if (keep.extra.columns & num_extra_columns > 0) {
    print(paste("Keeping", num_extra_columns, "extra columns"))
    colnames(output) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attribute",
                          paste0("extra_col", c(1:num_extra_columns)))
  } else if (num_extra_columns > 0) {
    print(paste("Removing", num_extra_columns, "extra columns"))
    output <- output[,c(1:9)]
    colnames(output) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attribute")
  } else {
    colnames(output) <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "attribute")
  }
  output # this returns the gtf as a data.frame
}



####
## Creating gene, transcript, and exon annotation files

# get_attr_list converts the attribute field of a GTF into a list, one element per row of GTF
get_attr_list <- function(gtf.df) {
  # remove unnecessary characters
  full_attr <- gsub("; ", ";", as.character(gtf.df$attribute))
  full_attr <- gsub("\"", "", full_attr)
  # split into a list of attributes
  full_attr_list <- strsplit(full_attr, ";")
  full_attr_list
}


# pull_attr gets one particular attribute from a list of attributes (list length = nrow(gtf), 
# each element of the list is a character vector of all of that row's attributes; 
# e.g. "gene_id ENSG00000223972.5" "transcript_id ENST00000450305.2" "gene_type transcribed_unprocessed_pseudogene")
pull_attr <- function(attr_list, pull_term) {
  outvals <- lapply(attr_list, grep, pattern = pull_term, value = T)
  outvals <- unlist(outvals)
  outvals <- gsub(paste0(pull_term, " "), "", outvals)
  outvals
}


# get_basic_attr pulls the basic gene info from the attributes field of a gtf
# input is a gtf; the attribute field is then prepped into a list that will be the
# input for pull_attr()
get_basic_attr <- function(gtf.df) {
  full_attr_list <- get_attr_list(gtf.df)
  output <- cbind(gtf.df, ensgene = pull_attr(full_attr_list, "gene_id"), 
                  symbol = pull_attr(full_attr_list, "gene_name"), biotype = pull_attr(full_attr_list, "gene_type"), 
                  gencode_level = paste0("level_", pull_attr(full_attr_list, "^level")))
  output
}


# get_tx_attr pulls additional transcript-specific info from a gtf file
# (can be used on top of get_basic_attr())
get_tx_attr <- function(gtf.df) {
  output <- gtf.df[gtf.df$type == "transcript", ]
  full_attr_list <- get_attr_list(output)
  output <- cbind(output, tx_id = pull_attr(full_attr_list, "transcript_id"),
                  tx_type = pull_attr(full_attr_list, "transcript_type"))
  output
}


# get_exon_attr pulls additional transcript and exon-specific info from a gtf file
# (can be used on top of get_basic_attr())
get_exon_attr <- function(gtf.df) {
  output <- gtf.df[gtf.df$type == "exon", ]
  full_attr_list <- get_attr_list(output)
  output <- cbind(output, tx_id = pull_attr(full_attr_list, "transcript_id"),
                  tx_type = pull_attr(full_attr_list, "transcript_type"), 
                  exon_number = as.numeric(pull_attr(full_attr_list, "exon_number")))
  output
}

