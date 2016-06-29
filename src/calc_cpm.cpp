#include "scater.h"

template <typename T>
SEXP calc_cpm_internal (const T* ptr, const matrix_info& MAT, SEXP lib_sizes, SEXP prior_counts, SEXP log) {
    if (!isReal(lib_sizes) || LENGTH(lib_sizes)!=MAT.ncol) { 
        throw std::runtime_error("length of 'lib_sizes' does not equal number of columns");
    }
    const double* lptr=REAL(lib_sizes);
    if (!isReal(prior_counts) || LENGTH(prior_counts)!=MAT.ncol) { 
        throw std::runtime_error("length of 'prior_counts' does not equal number of columns");
    }
    const double* pptr=REAL(prior_counts);
    if (!isLogical(log) || LENGTH(log)!=1) {
        throw std::runtime_error("log specification should be a logical scalar"); 
    }
    const bool dolog=asLogical(log);

    SEXP output=PROTECT(allocMatrix(REALSXP, MAT.nrow, MAT.ncol));
    try {
        const T* original=ptr;
        double* optr=REAL(output);
        for (size_t c=0; c<MAT.ncol; ++c) {
            for (size_t r=0; r<MAT.nrow; ++r) {
                optr[r]= (ptr[r] + pptr[c])/lptr[c];
            }
            ptr+=MAT.nrow;
            optr+=MAT.nrow;
        }

        if (dolog) {
            ptr=original;
            optr=REAL(output);
            for (size_t c=0; c<MAT.ncol; ++c) {
                for (size_t r=0; r<MAT.nrow; ++r) {
                    optr[r]=std::log(optr[r])/M_LN2;
                }
                ptr+=MAT.nrow;
                optr+=MAT.nrow;
            }
        }
    } catch (std::exception &e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}

SEXP calc_cpm(SEXP counts, SEXP lib_sizes, SEXP prior_counts, SEXP log) try {
    matrix_info MAT=check_matrix(counts);
    if (MAT.is_integer) {
        return calc_cpm_internal<int>(MAT.iptr, MAT, lib_sizes, prior_counts, log);
    } else {
        return calc_cpm_internal<double>(MAT.dptr, MAT, lib_sizes, prior_counts, log);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}
