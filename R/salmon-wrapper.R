## Wrappers for Salmon quantification of abundance of RNA-seq reads


################################################################################
#' Salmon wrapper functions
#' 
#' @param directory character string giving the path to the directory containing
#' the Salmon results for the sample. 
#' 
#' @details The directory is expected to contain results for just a single 
#' sample. Putting more than one sample's results in the directory will result
#' in unpredictable behaviour with this function. The function looks for the 
#' files (with the default names given by Salmon) 'quant.sf', 
#' 'stats.tsv', 'libFormatCounts.txt' and the sub-directories 'logs' (which 
#' contains a log file) and 'libParams' (which contains a file detailing the 
#' fragment length distribution). If these files are missing, or if results 
#' files have different names, then this function will not find them. 
#' 
#' This function will work for Salmon v0.7.x and greater, as the name of one of
#' the default output directories was changed from "aux" to "aux_info" in 
#' Salmon v0.7.
#' 
#' @return A list with two elements: (1) a data.frame \code{abundance} with 
#' columns for 'target_id' (feature, transcript, gene etc), 'length' (feature 
#' length), 'est_counts' (estimated feature counts), 'tpm' (transcripts per 
#' million); (2) a list, \code{run_info}, with metadata about the Salmon run that 
#' generated the results, including number of reads processed, mapping 
#' percentage, the library type used for the RNA-sequencing, including details about 
#' number of reads that did not match the given or inferred library type,
#' details about the Salmon command used to generate the results, and so on.
#' 
#' @name salmon-wrapper
#' @rdname salmon-wrapper
#' @export
#' @examples
#' \dontrun{
#' # If Salmon results are in the directory "output", then call:
#' readSalmonResultsOneSample("output")
#' }
readSalmonResultsOneSample <- function(directory) {
    cat(directory, "\n")
    ## Read in abundance information for the sample
    file_to_read <- file.path(directory, "quant.sf")
    abundance <- run_info <- NULL
    if ( file.exists(file_to_read) ) {
        ## read abundance values
        abundance <- data.table::fread(file_to_read, sep = "\t")
        abundance <- as.data.frame(abundance, stringsAsFactors = FALSE)
        colnames(abundance) <- c("target_id", "length", "eff_length", "tpm", "est_counts")
        
        ## extract run info
        json_file <- file.path(directory, "aux_info", "meta_info.json")
        jnames <- c("salmon_version", "samp_type", "num_libraries",
                    "library_types", "frag_dist_length", "seq_bias_correct",
                    "gc_bias_correct", "num_bias_bins", "mapping_type",
                    "num_targets", "num_bootstraps", "num_processed",
                    "num_mapped", "percent_mapped", "start_time")
        if (!file.exists(json_file))
            stop(paste(json_file, "not found or does not exist."))
        run_info <- as.list(rep(NA, length(jnames)))
        names(run_info) <- jnames
                tryCatch({
            tmp <- rjson::fromJSON(file = json_file)
            for (x in jnames) {
                if (x %in% names(tmp))
                    run_info[[x]] <- tmp[[x]]
            }
        }, error = function(e) {
            cat(paste("Reading aux_info/meta_info.json failed for", directory, "\n"))
        }, finally = {})
        
        ## extract library format info
        libformat_json_file <- file.path(directory, "lib_format_counts.json")
        if (!file.exists(libformat_json_file)) 
            stop(paste(libformat_json_file, "not found or does not exist."))
        lnames <- c("read_files", "expected_format", "compatible_fragment_ratio", 
                    "num_compatible_fragments", "num_assigned_fragments", 
                    "num_consistent_mappings", "num_inconsistent_mappings",
                    "strand_mapping_bias")
        libformat_info <- as.list(rep(NA, length(lnames)))
        names(libformat_info) <- lnames
        tryCatch({
            tmp <- rjson::fromJSON(file = libformat_json_file)
            for (x in lnames) {
                if (x %in% names(tmp))
                    libformat_info[[x]] <- tmp[[x]]
            }
        }, error = function(e) {
            cat(paste("Reading lib_format_counts.json failed for", directory, "\n"))
        }, finally = {})
        
        ## extract command info
        cmd_json_file <- file.path(directory, "cmd_info.json")
        if (!file.exists(cmd_json_file)) 
            stop(paste(cmd_json_file, "not found or does not exist."))
        cnames <- c("index", "libType", "mates1", "mates2", "output", "auxDir")
        cmd_info <- as.list(rep(NA, length(cnames)))
        names(cmd_info) <- cnames
        tryCatch({
            cmd_info <- rjson::fromJSON(file = cmd_json_file)[cnames]
        }, error = function(e) {
            cat(paste("Reading cmd_info.json failed for", directory, "\n"))
        }, finally = {})

        ## combine run info
        info_out <- c(run_info, libformat_info, cmd_info)
        
    } else
        stop(paste("File", file_to_read, "not found or does not exist. Please check directory is correct."))
    ## output list with abundances and run info
    list(abundance = abundance, run_info = info_out)
}


################################################################################
#' Read Salmon results from a batch of jobs 
#' 
#' After generating transcript/feature abundance results using Salmon for a 
#' batch of samples, read these abundance values into a
#'  \code{SingleCellExperiment} object.
#' 
#' @param Salmon_log list, generated by \code{runSalmon}. If provided, then 
#' \code{samples} and \code{directories} arguments are ignored.
#' @param samples character vector providing a set of sample names to use for 
#' the abundance results.
#' @param directories character vector providing a set of directories containing
#' Salmon abundance results to be read in.
#' @param logExprsOffset numeric scalar, providing the offset used when doing
#' log2-transformations of expression data to avoid trying to take logs of zero.
#' Default offset value is \code{1}.
#' @param verbose logical, should function provide output about progress?
#' 
#' @details This function expects to find only one set of Salmon abundance 
#' results per directory; multiple adundance results in a given directory will 
#' be problematic.
#' 
#' @return an SingleCellExperiment object
#' 
#' @name salmon-wrapper
#' @rdname salmon-wrapper
#' @export
#' @examples
#' \dontrun{
#' ## Define output directories in a vector called here "Salmon_dirs"
#' ## and sample names as "Salmon_samples"
#' sceset <- readSalmonResults(samples = Salmon_samples, 
#' directories = Salmon_dirs)
#' }
#' 
readSalmonResults <- function(Salmon_log = NULL, samples = NULL, 
                              directories = NULL, logExprsOffset = 1, 
                              verbose = TRUE) {
    ## initialise failure vector
    Salmon_fail <- rep(FALSE, length(samples))
    ## Checks on arguments
    if ( !is.null(Salmon_log) ) {
        cat("Using Salmon_log argument to define samples and results directories.")
        if ( !is.list(Salmon_log) )
            stop("The Salmon_log argument should be a list returned by runSalmon()")
        samples <- names(Salmon_log)       
        directories <- sapply(Salmon_log, function(x) {x$output_dir})
        logs <- lapply(Salmon_log, function(x) {x$Salmon_log})
        ## Can only check Salmon fail if log provided
        Salmon_fail <- sapply(logs, function(x) {
            any(grepl("[wW]arning|[eE]rror", x))})
        if ( any(Salmon_fail) ) {
            warning(paste0("The Salmon job failed for the following samples:\n ",
                           paste0(names(logs)[Salmon_fail], collapse = "\n"),
                           "\n It is recommended that you inspect Salmon_log for these samples."))
        }
        
    } else {
        cat("Salmon log not provided - assuming all runs successful")
        if ( is.null(samples) | is.null(directories) )
            stop("If Salmon_log argument is not used, then both samples and directories must be provided.")
        if ( length(samples) != length(directories) )
            stop("samples and directories arguments must be the same length")
    }
    
    samples <- samples[!Salmon_fail]
    directories <- directories[!Salmon_fail]
    
    ## Read first file to get size of feature set
    s1 <- readSalmonResultsOneSample(directories[1])
    nsamples <- length(samples)
    nfeatures <- nrow(s1$abundance)
    info_vars <- names(s1$run_info)
    ninfo_vars <- length(info_vars)
        
    ## Currently not reading in bootstrap results - to add in future
    
    ## Set up results objects
    pdata <- data.frame(matrix(NA, nrow = nsamples, ncol = ninfo_vars))
    rownames(pdata) <- samples
    colnames(pdata) <- names(s1$run_info)
    fdata <- data.frame(feature_id = s1$abundance$target_id,
                        feature_length = s1$abundance$length,
                        stringsAsFactors = FALSE)
    rownames(fdata) <- s1$abundance$target_id
    est_counts <- tpm <- feat_eff_len <- 
        matrix(NA, nrow = nfeatures, ncol = nsamples)
    colnames(est_counts) <- colnames(tpm) <- colnames(feat_eff_len) <- samples
    rownames(est_counts) <- rownames(tpm) <- rownames(feat_eff_len) <- 
        s1$abundance$target_id
    
    ## Read Salmon results into results objects
    if ( verbose )
        cat(paste("\nReading results for", nsamples, "samples:\n"))
    for (i in seq_len(nsamples)) {
        tmp_samp <- readSalmonResultsOneSample(directories[i])
        ## counts
        if ( length(tmp_samp$abundance$est_counts) != nfeatures )
            warning(paste("Results for directory", directories[i], 
                          "do not match dimensions of other samples."))
        else            
            est_counts[, i] <- tmp_samp$abundance$est_counts
        ## tpm
        if ( length(tmp_samp$abundance$est_counts) == nfeatures )
            tpm[, i] <- tmp_samp$abundance$tpm
        ## feature effective length
        if ( length(tmp_samp$abundance$eff_length) == nfeatures )
            feat_eff_len[, i] <- tmp_samp$abundance$eff_length
        ## run info
        for (x in info_vars)
            if (x %in% names(tmp_samp$run_info))
                pdata[[x]][i] <- tmp_samp$run_info[[x]]
        ## in future, add code to read in bootstraps
        if ( verbose ) {
            cat(".")
            if ( i %% 80 == 0)
                cat("\n")
        }
    }
    ## Add median feature effective length to fData
    fdata$median_effective_length <- matrixStats::rowMedians(feat_eff_len)
    
    if ( verbose )
        cat("\n")
    ## Produce SCESet object
    pdata <- pdata[!duplicated(rownames(pdata)), !duplicated(colnames(pdata))]
    fdata <- fdata[!duplicated(rownames(fdata)), !duplicated(colnames(fdata))]
    ## Produce SingleCellExperiment object
    sce_out <- SingleCellExperiment(
        list(exprs = log2(tpm + logExprsOffset), 
             counts = est_counts, tpm = tpm, 
             feature_effective_length = feat_eff_len), 
        colData = pdata, rowData = fdata)
    if ( verbose )
        cat("Using log2(TPM + 1) as 'exprs' values in output.\n")
    ## Return SCESet object
    sce_out
}


################################################################################
#' Run Salmon on FASTQ files to quantify feature abundance
#'
#' Run the abundance quantification tool \code{Salmon} on a set of FASTQ
#' files. Requires \code{Salmon} (\url{https://combine-lab.github.io/salmon/})
#' to be installed and a Salmon transcript index must have been generated prior
#' to using this function. See the Salmon website for installation and basic
#' usage instructions.
#'
#' @param targets_file character string giving the path to a tab-delimited text
#' file with either 2 columns (single-end reads) or 3 columns (paired-end reads)
#' that gives the sample names (first column) and FastQ file names (column 2 and
#' if applicable 3). The file is assumed to have column headers, although these
#' are not used.
#' @param transcript_index character string giving the path to the Salmon
#' index to be used for the feature abundance quantification.
#' @param single_end logical, are single-end reads used, or paired-end reads?
#' @param output_prefix character string giving the prefix for the output folder
#' that will contain the Salmon results. The default is \code{"output"} and
#' the sample name (column 1 of \code{targets_file}) is appended (preceded by an
#' underscore).
#' @param lib_type scalar, indicating RNA-seq library type. See Salmon 
#' documentation for details. Default is "A", for automatic detection.
#' @param n_processes integer giving the number of processes to use for
#' parallel Salmon jobs across samples. The package \code{parallel} is used. 
#' Default is 2 concurrent processes.
#' @param n_thread_per_process integer giving the number of threads for Salmon
#' to use per process (to parallelize Salmon for a given sample). Default is 4.
#' @param n_bootstrap_samples integer giving the number of bootstrap samples
#' that Salmon should use (default is 0). With bootstrap samples, uncertainty
#' in abundance can be quantified.
#' @param seqBias logical, should Salmon's option be used to model and correct
#' abundances for sequence specific bias? Default is \code{TRUE}.
#' @param gcBias logical, should Salmon's option be used to model and correct
#' abundances for GC content bias? Requires Salmon version 0.7.2 or higher. 
#' Default is \code{TRUE}.
#' @param posBias logical, should Salmon's option be used to model and correct
#' abundances for positional biases? Requires Salmon version 0.7.3 or higher.
#' Default is \code{FALSE}.
#' @param allowOrphans logical, Consider orphaned reads as valid hits when 
#' performing lightweight-alignment. This option will increase sensitivity 
#' (allow more reads to map and more transcripts to be detected), but may 
#' decrease specificity as orphaned alignments are more likely to be spurious.
#' For more details see Salmon documentation.
#' @param advanced_opts character scalar supplying list of advanced option 
#' arguments to apply to each Salmon call. For details see Salmon documentation
#' or type \code{salmon quant --help-reads} at the command line.
#' @param dry_run logical, if \code{TRUE} then a list containing the Salmon
#' commands that would be run and the output directories is returned. Can be
#' used to read in results if Salmon is run outside an R session or to produce
#'  a script to run outside of an R session.
#' @param salmon_cmd (optional) string giving full command to use to call
#' Salmon, if simply typing "salmon" at the command line does not give the
#' required version of Salmon or does not work. Default is simply "salmon".
#' If used, this argument should give the full path to the desired Salmon
#' binary.
#'
#' @details A Salmon transcript index can be built from a FASTA file:
#' \code{salmon index [arguments] FASTA-file}. See the Salmon documentation
#' for further details. This simple wrapper does not give access to all nuances
#' of Salmon usage. For finer-grained usage of Salmon please run it at the 
#' command line - results can still be read into R with 
#' \code{\link{readSalmonResults}}.
#'
#' @return A list containing three elements for each sample for which feature
#' abundance has been quantified: (1) \code{salmon_call}, the call used for
#' Salmon, (2) \code{salmon_log} the log generated by Salmon, and (3)
#' \code{output_dir} the directory in which the Salmon results can be found.
#'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapply
#' @importFrom utils read.delim
#' 
#' @name salmon-wrapper
#' @rdname salmon-wrapper
#' @export
#' @examples
#' \dontrun{
#' ## If in Salmon's 'test' directory, then try these calls:
#' ## Generate 'targets.txt' file:
#' write.table(data.frame(Sample="sample1", File1="reads_1.fastq.gz", File2="reads_1.fastq.gz"),
#'  file="targets.txt", quote=FALSE, row.names=FALSE, sep="\t")
#' Salmon_log <- runSalmon("targets.txt", "transcripts.idx", single_end=FALSE,
#'          output_prefix="output", verbose=TRUE, n_bootstrap_samples=10,
#'          dry_run = FALSE)
#' }
runSalmon <- function(targets_file, transcript_index, single_end = FALSE,
                      output_prefix = "output", lib_type = "A",
                      n_processes = 2, n_thread_per_process = 4,
                      n_bootstrap_samples = 0,
                      seqBias = TRUE, gcBias = TRUE, posBias = FALSE, 
                      allowOrphans = FALSE,
                      advanced_opts = NULL,
                      verbose = TRUE, dry_run = FALSE,
                      salmon_cmd = "salmon") {
    targets_dir <- paste0(dirname(targets_file), "/")
    targets <- read.delim(targets_file, stringsAsFactors = FALSE, header = TRUE)
    if ( !(ncol(targets) == 2 || ncol(targets) == 3) )
        stop("Targets file must have either 2 columns (single-end reads) or
             3 columns (paired-end reads). File should be tab-delimited with
             column headers")
    if ( ncol(targets) == 2 && !single_end ) {
        warning("targets only has two columns; proceeding assuming single-end
                reads")
        single_end <- TRUE
    }
    if ( ncol(targets) == 3 && single_end ) {
        warning("targets only has three columns, but 'single_end' was TRUE;
                proceeding assuming paired-end reads")
        single_end <- FALSE
    }
    samples <- targets[,1]
    ## Make sure that we'll be able to find the fastq files
    targets[, 2] <- paste0(targets_dir, targets[, 2])
    if ( ncol(targets) == 3 )
        targets[, 3] <- paste0(targets_dir, targets[, 3])
    ## Generate calls to Salmon
    output_dirs <- paste(output_prefix, samples, sep = "_")
    salmon_args <- paste("quant -i", transcript_index, "-o", output_dirs, 
                         "-l", lib_type, 
                         "--threads", n_thread_per_process,
                         "--numBootstraps", n_bootstrap_samples)
    names(salmon_args) <- samples
    if ( seqBias )
        salmon_args <- paste0(salmon_args, " --seqBias")
    if ( gcBias )
        salmon_args <- paste0(salmon_args, " --gcBias")
    if ( posBias )
        salmon_args <- paste0(salmon_args, " --posBias")
    if ( single_end )
        salmon_args <- paste(salmon_args, "-r", targets[, 2])
    else
        salmon_args <- paste(salmon_args, "-1", targets[, 2], 
                             "-2", targets[, 3])
    if ( allowOrphans )
        salmon_args <- paste(salmon_args, "--allowOrphans")
    if ( !is.null(advanced_opts) )
        salmon_args <- paste(salmon_args, advanced_opts)
    ##
    if ( dry_run ) {
        salmon_log <- vector("list", length(samples))
        names(salmon_log) <- samples
        Salmon_calls <- paste("salmon", salmon_args)
        for (i in seq_len(length(Salmon_calls))) {
            salmon_log[[i]]$output_dir <- output_dirs[i]
            salmon_log[[i]]$Salmon_call <- Salmon_calls[i]
            salmon_log[[i]]$salmon_log <-
                "Dry run: Salmon commands not executed."
        }
    } else {
        if (verbose)
            print(paste("Analysis started: ", Sys.time()))
        cl <- parallel::makeCluster(n_processes)
        # one or more parLapply calls to Salmon
        salmon_log <- parallel::parLapply(cl, salmon_args, .call_Salmon,
                                          salmon_cmd, verbose)
        parallel::stopCluster(cl)
        ## Return log of Salmon jobs, so user knows where to find results
        names(salmon_log) <- samples
        if (verbose) {
            print(paste("Analysis completed: ", Sys.time()))
            print(paste("Processed", length(samples), "samples"))
        }
        for (i in seq_len(length(salmon_log))) {
            salmon_log[[i]]$output_dir <- output_dirs[i]
        }
    }
    salmon_log
}

.call_Salmon <- function(scall, salmon_cmd, verbose = TRUE) {
    out <- tryCatch(ex <- system2(salmon_cmd, scall, stdout = TRUE,
                                  stderr = TRUE),
                    warning = function(w){w}, error = function(e){e})
    list(Salmon_call = paste(salmon_cmd, scall), salmon_log = out)
}




