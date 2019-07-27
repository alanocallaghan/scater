# Internal utilies related to the processing of size factors.

.get_all_sf_sets <- function(object) 
## Returns a list containing a list of sets of size factor values,
## as well as the indices of the rows to which each size factor set is to be applied.
{
    fcontrols <- spikeNames(object)
    
    # Storing the default size factors.
    sf.list <- vector("list", length(fcontrols)+1)
    sf.list[[1]] <- sizeFactors(object)
    to.use <- rep(1L, nrow(object))

    # Filling up the controls.
    counter <- 1L
    for (fc in fcontrols) {
        specific_sf <- sizeFactors(object, type=fc)
        if (is.null(specific_sf)) {
            warning(sprintf("spike-in set '%s' should have its own size factors", fc))
        } else {
            counter <- counter+1L
            which.current <- isSpike(object, type=fc)
            to.use[which.current] <- counter # after increment, as 1 is the NULL sizeFactors.
            sf.list[[counter]] <- specific_sf
        }
    }

    # Returning the output.
    return(list(size_factors=sf.list[seq_len(counter)], index=to.use)) 
}

.apply_to_size_factors <- function(object, FUN) 
## Running through the sets of size factors and centering them as necessary.
{
    sf <- sizeFactors(object)
    if (!is.null(sf)) {
        sizeFactors(object) <- FUN(sf)
    }
     
    for (sf_name in sizeFactorNames(object)) {
        sf <- sizeFactors(object, sf_name)
        if (!is.null(sf)) {
            sizeFactors(object, sf_name) <- FUN(sf)
        }
    }

    return(object)
}

.replace_size_factors <- function(object, use_size_factors) 
## Either eliminates all size factors in a SingleCellExperiment object,
## or sets the main size factors to the specified vector.
{
    if (is.logical(use_size_factors)) {
        if (!use_size_factors) {
            object <- clearSizeFactors(object)

            # Also eliminating spike-in information, to avoid warnings
            # when spike-in size factors are not available.
            object <- clearSpikes(object)
        }
    } else {
        sizeFactors(object) <- rep(use_size_factors, length.out=ncol(object))
    }
    return(object)
}
