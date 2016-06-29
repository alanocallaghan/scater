#include "scater.h"

bool isNA(int x) {
    return x==NA_INTEGER;
}

bool isNA(double x) {
    return ISNA(x);
}

/*********************************************************
 This contains a number of small functions, written to improve 
 speed or memory efficiency over a native R implementation. 
 **********************************************************/

/* A function to get the column sum in a subset of rows. */

template <typename T>
SEXP colsum_subset_internal (const T* ptr, const matrix_info& MAT, SEXP subset) {
    if (!isInteger(subset)) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }
    const int slen=LENGTH(subset);
    const int* sptr=INTEGER(subset);
    for (int s=0; s<slen; ++s) {
        if (sptr[s] < 1 || sptr[s] > MAT.nrow) { 
            throw std::runtime_error("subset indices out of range");
        }
    }

    SEXP output=PROTECT(allocVector(REALSXP, MAT.ncol));
    try {
        double* optr=REAL(output);
        
        // Summing across, using 1-indexed pointers.
        --ptr;
        int s;
        for (size_t c=0; c<MAT.ncol; ++c) {
            optr[c]=0;
            for (s=0; s<slen; ++s) {
                optr[c]+=ptr[sptr[s]];
            }
            ptr+=MAT.nrow;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);
    return output;
}

SEXP colsum_subset(SEXP matrix, SEXP subset) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return colsum_subset_internal<int>(MAT.iptr, MAT, subset);
    } else {
        return colsum_subset_internal<double>(MAT.dptr, MAT, subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* A function to get the number of above-threshold values in each column, for a subset of rows. */

template <typename T>
SEXP colsum_exprs_subset_internal (const T* ptr, const matrix_info& MAT, T threshold, SEXP subset) {
    if (!isInteger(subset)) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }
    const int slen=LENGTH(subset);
    const int* sptr=INTEGER(subset);
    for (int s=0; s<slen; ++s) {
        if (sptr[s] < 1 || sptr[s] > MAT.nrow) { 
            throw std::runtime_error("subset indices out of range");
        }
    }

    SEXP output=PROTECT(allocVector(INTSXP, MAT.ncol));
    try {
        int* optr=INTEGER(output);
        
        // Summing across, using 1-indexed pointers.
        --ptr;
        int s;
        for (size_t c=0; c<MAT.ncol; ++c) {
            optr[c]=0;
            for (s=0; s<slen; ++s) {
                if (ptr[sptr[s]] > threshold) {
                    ++(optr[c]);
                }
            }
            ptr+=MAT.nrow;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);
    return output;
}

SEXP colsum_exprs_subset(SEXP matrix, SEXP threshold, SEXP subset) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        if (!isInteger(threshold) || LENGTH(threshold)!=1) { 
            throw std::runtime_error("threshold should be an integer scalar");
        }
        return colsum_exprs_subset_internal<int>(MAT.iptr, MAT, asInteger(threshold), subset);
    } else {
        if (!isReal(threshold) || LENGTH(threshold)!=1) { 
            throw std::runtime_error("threshold should be a double-precision scalar");
        }
        return colsum_exprs_subset_internal<double>(MAT.dptr, MAT, asReal(threshold), subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* A function to get the number of above-threshold values in each row. */

template <typename T>
SEXP rowsum_exprs_internal (const T* ptr, const matrix_info& MAT, T threshold) {
    SEXP output=PROTECT(allocVector(INTSXP, MAT.nrow));
    try {
        int* optr=INTEGER(output);
        for (size_t r=0; r<MAT.nrow; ++r) {
            optr[r]=0;
        }
        
        // Summing across, using 1-indexed pointers.
        for (size_t c=0; c<MAT.ncol; ++c) {
            for (size_t r=0; r<MAT.nrow; ++r) {
                if (ptr[r] > threshold) {
                    ++(optr[r]);
                }
            }
            ptr+=MAT.nrow;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }
    UNPROTECT(1);
    return output;
}

SEXP rowsum_exprs(SEXP matrix, SEXP threshold) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        if (!isInteger(threshold) || LENGTH(threshold)!=1) { 
            throw std::runtime_error("threshold should be an integer scalar");
        }
        return rowsum_exprs_internal<int>(MAT.iptr, MAT, asInteger(threshold));
    } else {
        if (!isReal(threshold) || LENGTH(threshold)!=1) { 
            throw std::runtime_error("threshold should be a double-precision scalar");
        }
        return rowsum_exprs_internal<double>(MAT.dptr, MAT, asReal(threshold));
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* A function to check whether the counts are non-negative or missing. */

template <typename T>
SEXP negative_counts_internal(const T* ptr, const matrix_info& MAT) {
    const size_t total_size=MAT.nrow*MAT.ncol;
    for (size_t i=0; i<total_size; ++i) {
        if (ptr[i] < 0 || isNA(ptr[i])) { return ScalarLogical(1); }
    } 
    return ScalarLogical(0);        
}

SEXP negative_counts(SEXP matrix) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return negative_counts_internal<int>(MAT.iptr, MAT);
    } else {
        return negative_counts_internal<double>(MAT.dptr, MAT);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* A function to check whether the expression values are NA. */

template <typename T>
SEXP missing_exprs_internal(const T* ptr, const matrix_info& MAT) {
    const size_t total_size=MAT.nrow*MAT.ncol;
    for (size_t i=0; i<total_size; ++i) {
        if (isNA(ptr[i])) { return ScalarLogical(1); }
    } 
    return ScalarLogical(0);        
}

SEXP missing_exprs(SEXP matrix) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return missing_exprs_internal<int>(MAT.iptr, MAT);
    } else {
        return missing_exprs_internal<double>(MAT.dptr, MAT);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

/* A function to get the percentage of counts/expression taken up by the top set of genes. */

template <typename T>
SEXP calc_top_features_internal(const T* ptr, const matrix_info& MAT, SEXP top) {
    if (!isInteger(top)) { 
        throw std::runtime_error("top specification must be an integer vector");
    }
    const int ntop=LENGTH(top);
    const int *tptr=INTEGER(top);
    for (size_t t=1; t<ntop; ++t) {
        if (tptr[t] < tptr[t-1]) { 
            throw std::runtime_error("numbers of top genes must be sorted"); 
        }
    }
    if (ntop && (tptr[0] < 1 || tptr[ntop-1] > MAT.nrow)) {
        throw std::runtime_error("number of top genes is out of index range");
    }
    
    SEXP output=PROTECT(allocMatrix(REALSXP, MAT.ncol, ntop));
    try {
        double** optrs=(double**)R_alloc(ntop, sizeof(double*));
        size_t t;
        if (ntop) {
            optrs[0]=REAL(output);
        }
        for (t=1; t<ntop; ++t) {
            optrs[t]=optrs[t-1]+MAT.ncol;
        }

        std::deque<T> values(MAT.nrow);
        size_t r, x;
        size_t target_index;
        double accumulated, total;

        for (size_t c=0; c<MAT.ncol; ++c) {
            total=0;
            for (r=0; r<MAT.nrow; ++r) {
                values[r]=ptr[r];
                total+=ptr[r];
            }
            
            // Sorting in descending order, and computing the accumulated total.
            std::sort(values.begin(), values.end(), std::greater<T>());
            x=0;
            accumulated=0;
            for (t=0; t<ntop; ++t) {
                target_index=size_t(tptr[t] - 1); // 0-based index.
                while (x<=target_index && x<MAT.nrow) {
                    accumulated+=values[x];
                    ++x;
                }
                optrs[t][c]=accumulated/total;
            }

            ptr+=MAT.nrow;
        }
    } catch (std::exception& e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP calc_top_features (SEXP matrix, SEXP top) try {
    matrix_info MAT=check_matrix(matrix);
    if (MAT.is_integer){
        return calc_top_features_internal<int>(MAT.iptr, MAT, top);
    } else {
        return calc_top_features_internal<double>(MAT.dptr, MAT, top);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}

