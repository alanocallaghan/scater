#include "scater.h"

template <typename T, class V, class M>
Rcpp::RObject calc_variance_internal(M mat, Rcpp::RObject subset, Rcpp::LogicalVector byrow) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    if (byrow.size()!=1) {
        throw std::runtime_error("'byrow' should be a logical scalar");
    }
    const bool BR=byrow[0];

    // Computing variances over rows or columns.
    if (BR) {
        V incoming(ncells);
        Rcpp::NumericVector outvar(ngenes);
        Rcpp::IntegerVector subout=process_subset_vector(subset, mat, false);

        if (!subout.size()) { 
            std::fill(outvar.begin(), outvar.end(), R_NaReal);
            return outvar;
        }

        for (size_t g=0; g<ngenes; ++g) {
            mat->get_row(g, incoming.begin());

            // Computing the mean.
            double curmean=0;
            for (const auto& s : subout) {
                curmean+=incoming[s];
            }
            curmean/=subout.size();

            // Computing the variance.
            double& curval=outvar[g];
            for (const auto& s : subout) { 
                const double tmp=incoming[s] - curmean;
                curval += tmp * tmp;
            }
            curval/=subout.size()-1;
        }

        return outvar;
    } else {
        V incoming(ngenes);
        Rcpp::NumericVector outvar(ncells);
        Rcpp::IntegerVector subout=process_subset_vector(subset, mat, true);

        if (!subout.size()) {
            std::fill(outvar.begin(), outvar.end(), R_NaReal);
            return outvar;
        }

        for (size_t c=0; c<ncells; ++c) {
            auto iIt=mat->get_const_col(c, incoming.begin());

            // Computing the mean.
            double curmean=0;
            for (const auto& s : subout) {
                curmean+=*(iIt+s);
            }
            curmean/=subout.size();

            // Computing the variance.
            double& curval=outvar[c];
            for (const auto& s : subout) { 
                const double tmp=*(iIt+s) - curmean;
                curval += tmp*tmp;
            }
            curval/=subout.size()-1;
        }

        return outvar;
    }
}

SEXP calc_variance(SEXP exprs, SEXP subset, SEXP byrow) {
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return calc_variance_internal<int, Rcpp::IntegerVector>(mat.get(), subset, byrow);

    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return calc_variance_internal<double, Rcpp::NumericVector>(mat.get(), subset, byrow);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

