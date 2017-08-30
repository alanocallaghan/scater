#include "scater.h"

template <typename T, class V, class M>
Rcpp::RObject margin_summary_internal(M mat, T threshold, Rcpp::RObject subset, Rcpp::LogicalVector byrow) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    if (byrow.size()!=1) {
        throw std::runtime_error("'byrow' should be a logical scalar");
    }
    const bool BR=byrow[0];

    // Computing margin summaries over rows or over columns (using column extraction in both cases).
    V incoming(ngenes);
    if (BR) {
        V outsum(ngenes);
        Rcpp::IntegerVector outnum(ngenes);
        Rcpp::IntegerVector subout=process_subset_vector(subset, mat, false);

        for (const auto& s : subout) {
            auto iIt=mat->get_const_col(s, incoming.begin());
            
            auto osIt=outsum.begin();
            auto onIt=outnum.begin();
            for (size_t r=0; r<ngenes; ++r) {
                const auto& current=*(iIt+r);
                *osIt += current;
                if (current > threshold) { ++(*onIt); }
                ++osIt;
                ++onIt;
            }
        }

        return Rcpp::List::create(outsum, outnum);
    } else {
        V outsum(ncells);
        Rcpp::IntegerVector outnum(ncells);
        Rcpp::IntegerVector subout=process_subset_vector(subset, mat, true);

        auto osIt=outsum.begin();
        auto onIt=outnum.begin();
        for (size_t c=0; c<ncells; ++c) {
            auto iIt=mat->get_const_col(c, incoming.begin());

            auto& cursum=(*osIt);
            auto& curnum=(*onIt);
            for (const auto& s : subout) {

                const auto& current=*(iIt+s);
                cursum += current;
                if (current > threshold) { ++curnum; }
            }

            ++osIt;
            ++onIt;
        }

        return Rcpp::List::create(outsum, outnum);
    }
}

SEXP margin_summary (SEXP exprs, SEXP threshold, SEXP subset, SEXP byrow) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        Rcpp::IntegerVector thresh(threshold);
        if (thresh.size()!=1) {
            throw std::runtime_error("threshold should be an integer scalar");
        }
        return margin_summary_internal<int, Rcpp::IntegerVector>(mat.get(), thresh[0], subset, byrow);

    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        Rcpp::NumericVector thresh(threshold);
        if (thresh.size()!=1) {
            throw std::runtime_error("threshold should be a double-precision scalar");
        }
        return margin_summary_internal<double, Rcpp::NumericVector>(mat.get(), thresh[0], subset, byrow);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

