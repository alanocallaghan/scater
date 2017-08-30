#include "scater.h"

/* A function to get the percentage of counts/expression taken up by the top set of genes. */

template <typename T, class V, class M>
Rcpp::RObject calc_top_features_internal(M mat, Rcpp::RObject topvec, Rcpp::RObject subset) {

    if (topvec.sexp_type()!=INTSXP) {
        throw std::runtime_error("'top' must be an integer vector");
    }
    Rcpp::IntegerVector top(topvec);
    const size_t ntop=top.size();
    for (size_t t=1; t<ntop; ++t) {
        if (top[t] < top[t-1]) { 
            throw std::runtime_error("numbers of top genes must be sorted"); 
        }
    }

    // Checking subsetting vector.
    const bool use_subset=(subset!=R_NilValue);
    Rcpp::IntegerVector subvec;
    if (use_subset) {
        subvec=process_subset_vector(subset, mat, true);
    }

    const size_t used_genes=(use_subset ? subvec.size() : mat->get_nrow());
    if (ntop && (top[0] < 1 || top[ntop-1] > used_genes)) {
        throw std::runtime_error("number of top genes is out of index range");
    }
    
    // Setting up output vectors.
    const size_t ncells=mat->get_ncol(); 
    std::deque<Rcpp::NumericVector> temp_output(ntop);
    std::deque<Rcpp::NumericVector::iterator> oIts(ntop);
    for (size_t t=0; t<ntop; ++t) {
        temp_output[t]=Rcpp::NumericVector(ncells);
        oIts[t]=temp_output[t].begin();
    }
    V input(mat->get_nrow()), survivors(used_genes);

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, input.begin());

        // Subsetting.
        if (use_subset) {
            auto suIt=survivors.begin();
            for (auto svIt=subvec.begin(); svIt!=subvec.end(); ++svIt, ++suIt) { 
                (*suIt)=input[*svIt];
            }
        } else {
            survivors=input;
        }
        double total=std::accumulate(survivors.begin(), survivors.end(), 0);
            
        // Sorting in descending order, and computing the accumulated total.
        std::partial_sort(survivors.begin(), survivors.begin() + top[ntop-1], survivors.end(), std::greater<T>());
        size_t x=0;
        double accumulated=0;
        for (size_t t=0; t<ntop; ++t) { 
            size_t target_index=top[t] - 1; // 0-based index.
            while (x<=target_index && x<survivors.size()) {
                accumulated+=survivors[x];
                ++x;
            }
            *(oIts[t])=accumulated/total * 100;
            ++(oIts[t]);
        }
    }

    // Setting up output list.
    Rcpp::List output(ntop);
    for (size_t t=0; t<ntop; ++t) {
        output[t]=temp_output[t];
    }
    return output;
}

SEXP calc_top_features (SEXP matrix, SEXP top, SEXP subset) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(matrix);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(matrix);
        return calc_top_features_internal<int, Rcpp::IntegerVector>(mat.get(), top, subset);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(matrix);
        return calc_top_features_internal<double, Rcpp::NumericVector>(mat.get(), top, subset);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

