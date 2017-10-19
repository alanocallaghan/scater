#include "scater.h"

template <class M, class O>
void downsample_matrix_internal(M mat, O output, Rcpp::NumericVector prop) {
    const size_t ngenes=mat->get_nrow();
    Rcpp::IntegerVector incoming(ngenes), outgoing(ngenes);

    const size_t ncells=mat->get_ncol();
    if (prop.size()!=ncells) {
        throw std::runtime_error("length of 'prop' should be equal to number of cells");
    }
    Rcpp::RNGScope _rng;
    
    for (size_t i=0; i<ncells; ++i) {
        mat->get_col(i, incoming.begin());
        const double& curprop=prop[i];
        if (curprop < 0 || curprop > 1) { 
            throw std::runtime_error("downsampling proportion must lie in [0, 1]");
        }

        // Getting the total library size and the size to downsample.
        const int num_total=std::accumulate(incoming.begin(), incoming.end(), 0),
                  num_sample=std::round(curprop*num_total);

        // Setting up the output vector.
        std::fill(outgoing.begin(), outgoing.end(), 0);
        auto oIt=outgoing.begin();
        auto iIt=incoming.begin();
        int cumulative=0;
        if (incoming.size()) { 
            cumulative+=*iIt;
            ++iIt;
        }

        // Sampling scheme adapted from John D. Cook, https://stackoverflow.com/a/311716/15485.
        int current=0, num_selected=0;
        while (num_selected < num_sample) {
            const double u = unif_rand(); 

            if ( (num_total - current)*u < num_sample - num_selected) {
                // Current read is selected, we advance to that read's "index" (if we had instantiated the full [0, num_total) array).
                while (cumulative <= current) {
                    cumulative+=(*iIt);
                    ++iIt;
                    ++oIt;
                }
                ++(*oIt);
                ++num_selected;
            }
            ++current; 
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
        downsample_matrix_internal(mat.get(), out.get(), prop);
        return out->yield();
    } else {
        auto mat=beachmat::create_numeric_matrix(rmat);
        auto out=beachmat::create_numeric_output(mat->get_nrow(), mat->get_ncol(), otype);
        downsample_matrix_internal(mat.get(), out.get(), prop);
        return out->yield();
    }
    END_RCPP    
}
