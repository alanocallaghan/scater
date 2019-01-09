#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

#include <stdexcept>
#include <vector>

template <class V, class O, class M>
Rcpp::RObject sum_counts_internal(M mat, O out, const std::vector<Rcpp::IntegerVector>& feature_set) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();
    V holder_in(ngenes), holder_out(feature_set.size());

    for (size_t c=0; c<ncells; ++c) {
        auto cIt=mat->get_const_col(c, holder_in.begin());

        auto fsIt=feature_set.begin();
        for (auto oIt=holder_out.begin(); oIt!=holder_out.end(); ++oIt, ++fsIt) {
            const auto& curset=(*fsIt);
            auto& curval=(*oIt);

            curval=0;
            for (auto feat : curset) {
                curval+=*(cIt + feat);                
            }
        }
        out->set_col(c, holder_out.begin());
    }
    return out->yield();
}

SEXP sum_counts (SEXP counts, SEXP features) {
    BEGIN_RCPP
    Rcpp::List Features(features);
    std::vector<Rcpp::IntegerVector> feature_set(Features.size());
    const size_t nfeatures=feature_set.size();
    for (size_t fx=0; fx<nfeatures; ++fx) {
        feature_set[fx]=Features[fx];
    }

    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        const size_t ncells=mat->get_ncol();

        beachmat::output_param OPARAM(mat->get_matrix_type(), true, true);
        OPARAM.optimize_chunk_dims(nfeatures, ncells);
        auto out=beachmat::create_integer_output(nfeatures, ncells, OPARAM);

        return sum_counts_internal<Rcpp::IntegerVector>(mat.get(), out.get(), feature_set);

    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        const size_t ncells=mat->get_ncol();

        beachmat::output_param OPARAM(mat->get_matrix_type(), true, true);
        OPARAM.optimize_chunk_dims(nfeatures, ncells);
        auto out=beachmat::create_numeric_output(nfeatures, ncells, OPARAM);

        return sum_counts_internal<Rcpp::NumericVector>(mat.get(), out.get(), feature_set);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
