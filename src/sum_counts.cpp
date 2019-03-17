#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <vector>

template <class V, class O, class M>
Rcpp::RObject sum_counts_internal(M mat, O out, const std::vector<Rcpp::IntegerVector>& summable_set, size_t start_index, size_t end_index, bool sum_by_row) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    if (sum_by_row) {
        V holder_in(ngenes), holder_out(summable_set.size());
        if (end_index > ncells) {
            throw std::runtime_error("end index out of range");
        }

        for (size_t c=start_index; c<end_index; ++c) {
            mat->get_col(c, holder_in.begin());
    
            auto ssIt=summable_set.begin();
            for (auto oIt=holder_out.begin(); oIt!=holder_out.end(); ++oIt, ++ssIt) {
                const auto& curset=(*ssIt);
                auto& curval=(*oIt=0);
                for (auto feat : curset) {
                    curval+=holder_in[feat];
                }
            }
            out->set_col(c - start_index, holder_out.begin());
        }

    } else {
        V holder_in(ncells), holder_out(summable_set.size());
        if (end_index > ngenes) {
            throw std::runtime_error("end index out of range");
        }

        for (size_t g=start_index; g<end_index; ++g) {
            mat->get_row(g, holder_in.begin());
    
            auto ssIt=summable_set.begin();
            for (auto oIt=holder_out.begin(); oIt!=holder_out.end(); ++oIt, ++ssIt) {
                const auto& curset=(*ssIt);
                auto& curval=(*oIt=0);
                for (auto cell : curset) {
                    curval+=holder_in[cell];
                }
            }
            out->set_row(g - start_index, holder_out.begin());
        }
    }
    return out->yield();
}

SEXP sum_counts (SEXP counts, SEXP sumset, SEXP job_start, SEXP job_end, SEXP by_row) {
    BEGIN_RCPP
    Rcpp::List Sumset(sumset);
    std::vector<Rcpp::IntegerVector> summation_set(Sumset.size());
    const size_t nsummations=summation_set.size();
    for (size_t fx=0; fx<nsummations; ++fx) {
        summation_set[fx]=Sumset[fx];
    }

    // Determining the type and amount of work to do.
    bool sum_by_row=check_logical_scalar(by_row, "by-row specification");
    size_t start_index=check_integer_scalar(job_start, "start index");
    size_t end_index=check_integer_scalar(job_end, "end index");
    if (end_index < start_index) { throw std::runtime_error("start index is less than end index"); }
    
    // Determining the dimensionality of the output.
    const size_t niters=end_index - start_index;
    const size_t out_nrow=(sum_by_row ? nsummations : niters);
    const size_t out_ncol=(sum_by_row ? niters : nsummations);

    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        beachmat::output_param OPARAM(mat->get_class(), mat->get_package());
        auto out=beachmat::create_integer_output(out_nrow, out_ncol, OPARAM);
        return sum_counts_internal<Rcpp::IntegerVector>(mat.get(), out.get(), summation_set, start_index, end_index, sum_by_row);

    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        beachmat::output_param OPARAM(mat->get_class(), mat->get_package());
        auto out=beachmat::create_numeric_output(out_nrow, out_ncol, OPARAM);
        return sum_counts_internal<Rcpp::NumericVector>(mat.get(), out.get(), summation_set, start_index, end_index, sum_by_row);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
