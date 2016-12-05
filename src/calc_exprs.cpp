#include "scater.h"

template <typename T>
SEXP calc_exprs_internal (const T* ptr, const matrix_info& MAT, SEXP size_fac, SEXP prior_count, SEXP log, SEXP sum, SEXP subset) {
    if (!isReal(size_fac) || LENGTH(size_fac)!=MAT.ncol) { 
        throw std::runtime_error("length of 'size_fac' does not equal number of columns");
    }
    const double* szptr=REAL(size_fac);
    if (!isReal(prior_count) || LENGTH(prior_count)!=1) { 
        throw std::runtime_error("'prior_count' should be a logical scalar");
    }
    const double prior=asReal(prior_count);

    // Checking flags.
    if (!isLogical(log) || LENGTH(log)!=1) {
        throw std::runtime_error("log specification should be a logical scalar"); 
    }
    const bool dolog=asLogical(log);
    if (!isLogical(sum) || LENGTH(sum)!=1) {
        throw std::runtime_error("sum specification should be a sumical scalar"); 
    }
    const bool dosum=asLogical(sum);

    // Checking subset data
    subset_info subout=process_subset_vector(subset, MAT);
    const int slen=subout.first;
    const int* sptr=subout.second;

    SEXP output;
    if (dosum) { 
        output=PROTECT(allocVector(REALSXP, slen));
    } else {
        output=PROTECT(allocMatrix(REALSXP, slen, MAT.ncol));
    }
    try {
        double* optr=REAL(output);
        if (dosum) {
            std::fill(optr, optr+slen, 0);
        }
        double tmp;
        int s;

        for (size_t c=0; c<MAT.ncol; ++c) {
            for (s=0; s<slen; ++s) {
                tmp=ptr[sptr[s]]/szptr[c] + prior;
                if (dosum) { 
                    optr[s]+=tmp;
                } else if (dolog) { 
                    optr[s]=std::log(tmp)/M_LN2;
                } else {
                    optr[s]=tmp;
                }
            }
            ptr+=MAT.nrow;
            if (!dosum) { 
                optr+=slen;
            }
        }

        if (dosum && dolog) {
            for (s=0; s<slen; ++s) { 
                optr[s]=std::log(optr[s])/M_LN2;
            }
        }
    } catch (std::exception &e) {
        UNPROTECT(1);
        throw;
    }

    UNPROTECT(1);
    return output;
}


SEXP calc_exprs(SEXP counts, SEXP size_fac, SEXP prior_count, SEXP log, SEXP sum, SEXP subset) try {
    matrix_info MAT=check_matrix(counts);
    if (MAT.is_integer) {
        return calc_exprs_internal<int>(MAT.iptr, MAT, size_fac, prior_count, log, sum, subset);
    } else {
        return calc_exprs_internal<double>(MAT.dptr, MAT, size_fac, prior_count, log, sum, subset);
    }
} catch (std::exception& e) {
    return mkString(e.what());
}
