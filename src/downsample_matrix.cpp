#include "scater.h"

template <class V, class M, class O>
void downsample_matrix_internal(M mat, O output, Rcpp::NumericVector prop) {
    const size_t ngenes=mat->get_nrow();
    V incoming(ngenes), outgoing(ngenes);

    const size_t ncells=mat->get_ncol();
    if (prop.size()!=ncells) {
        throw std::runtime_error("length of 'prop' should be equal to number of cells");
    }
    Rcpp::RNGScope _rng;
    
    for (size_t i=0; i<ncells; ++i) {
        auto it=mat->get_const_col(i, incoming.begin());
        const double& curprop=prop[i];
        
        for (size_t j=0; j<ngenes; ++j) {
            const auto& val=*it;
            if (val > 0) { // i.e., non-zero.
                outgoing[j]=R::rbinom(val, curprop);
            } else {
                outgoing[j]=0;
            }
            ++it;
        }
        
        output->set_col(i, outgoing.begin());
    }

    /* Note that the RNGscope destructor may trigger a garbage collection.
     * I'm not sure that the object from output->yield() remains protected if a compiler does not implement RVO.
     * This would result in a copy and destruction, and a point in time at which the output memory is unprotected.
     */
    return;
}

SEXP downsample_matrix(SEXP rmat, SEXP prop) {
    BEGIN_RCPP
    int rtype=beachmat::find_sexp_type(rmat);
    auto otype=beachmat::output_param(rmat, false, true);
    if (rtype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(rmat);
        auto out=beachmat::create_integer_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal<Rcpp::IntegerVector>(mat.get(), out.get(), prop);
        return out->yield();
    } else {
        auto mat=beachmat::create_numeric_matrix(rmat);
        auto out=beachmat::create_numeric_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal<Rcpp::NumericVector>(mat.get(), out.get(), prop);
        return out->yield();
    }
    END_RCPP    
}
