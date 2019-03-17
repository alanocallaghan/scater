#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <vector>

std::vector<Rcpp::IntegerVector> create_summation(Rcpp::List Sumset) {
    const size_t nsummations=Sumset.size();
    std::vector<Rcpp::IntegerVector> summation_set(Sumset.size());
    for (size_t fx=0; fx<nsummations; ++fx) {
        summation_set[fx]=Sumset[fx];
    }
    return summation_set;
}

/***************************
 * Summing across features *
 ***************************/

template <class V, class O, class M>
Rcpp::RObject sum_row_counts_internal(M mat, O out, const std::vector<Rcpp::IntegerVector>& summable_set, size_t start_index, size_t end_index) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    V holder_in(ngenes), holder_out(summable_set.size());
    if (end_index > ncells) {
        throw std::runtime_error("end index out of range");
    }

    // Specializing for dense arrays to avoid copies.
    // Again, quite difficult to support sparsity due to the
    // need for arbitrary indexing of rows.
    auto raws=mat->set_up_raw();
    const bool is_dense=(mat->col_raw_type()=="dense");

    for (size_t c=start_index; c<end_index; ++c) {
        typename V::iterator it;
        if (is_dense) {
            mat->get_col_raw(c, raws);
            it=raws.get_values_start();

        } else {
            it=holder_in.begin();
            mat->get_col(c, it);
        }

        auto ssIt=summable_set.begin();
        for (auto& curval : holder_out) {
            const auto& curset=(*ssIt);
            curval=0;
            for (auto feat : curset) {
                curval += *(it + feat);
            }
            ++ssIt;
        }

        out->set_col(c - start_index, holder_out.begin());
    }
    return out->yield();
}

SEXP sum_row_counts (SEXP counts, SEXP sumset, SEXP job_start, SEXP job_end) {
    BEGIN_RCPP

    auto summation_set=create_summation(sumset);
    const size_t nsummations=summation_set.size();

    // Determining the type and amount of work to do.
    size_t start_index=check_integer_scalar(job_start, "start index");
    size_t end_index=check_integer_scalar(job_end, "end index");
    if (end_index < start_index) { throw std::runtime_error("start index is less than end index"); }
    const size_t ncells=end_index - start_index;

    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        beachmat::output_param OPARAM(mat->get_class(), mat->get_package());
        auto out=beachmat::create_integer_output(nsummations, ncells, OPARAM);
        return sum_row_counts_internal<Rcpp::IntegerVector>(mat.get(), out.get(), summation_set, start_index, end_index);

    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        beachmat::output_param OPARAM(mat->get_class(), mat->get_package());
        auto out=beachmat::create_numeric_output(nsummations, ncells, OPARAM);
        return sum_row_counts_internal<Rcpp::NumericVector>(mat.get(), out.get(), summation_set, start_index, end_index);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
} 

/************************
 * Summing across cells *
 ************************/

template <class V, class O, class M>
Rcpp::RObject sum_col_counts_internal(M mat, O out, const std::vector<Rcpp::IntegerVector>& summable_set, size_t start_index, size_t end_index) {
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    V holder_in(ngenes), holder_out(ngenes);
    if (end_index > ngenes) {
        throw std::runtime_error("end index out of range");
    }

    auto raws=mat->set_up_raw();
    const bool is_dense=(mat->col_raw_type()=="dense");
    const bool is_sparse=(mat->col_raw_type()=="sparse");

    size_t i_summed=0;
    for (const auto& curset : summable_set) {
        std::fill(holder_out.begin()+start_index, holder_out.begin()+end_index, 0);

        // Specialized for sparse matrices.
        if (is_sparse) {
            for (auto cell : curset) {
                mat->get_col_raw(cell, raws, start_index, end_index);
                auto n=raws.get_n();
                auto idx=raws.get_structure_start();
                auto val=raws.get_values_start();

                for (size_t i=0; i<raws.get_n(); ++i, ++idx, ++val) {
                    holder_out[*idx]+=*val;
                }
            }

        } else {
            for (auto cell : curset) {
                typename V::iterator it;
                if (is_dense) {
                    mat->get_col_raw(cell, raws);
                    it=raws.get_values_start();
                } else {
                    it=holder_in.begin();
                    mat->get_col(cell, it+start_index, start_index, end_index);
                }

                for (size_t i=start_index; i<end_index; ++i) {
                    holder_out[i] += *(it + i);
                }
            }
        }

        out->set_col(i_summed, holder_out.begin()+start_index);
        ++i_summed;
    }
    return out->yield();
}

SEXP sum_col_counts (SEXP counts, SEXP sumset, SEXP job_start, SEXP job_end) {
    BEGIN_RCPP

    auto summation_set=create_summation(sumset);
    const size_t nsummations=summation_set.size();

    // Determining the type and amount of work to do.
    size_t start_index=check_integer_scalar(job_start, "start index");
    size_t end_index=check_integer_scalar(job_end, "end index");
    if (end_index < start_index) { throw std::runtime_error("start index is less than end index"); }
    const size_t ngenes=end_index - start_index;

    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(counts);
        beachmat::output_param OPARAM(mat->get_class(), mat->get_package());
        auto out=beachmat::create_integer_output(ngenes, nsummations, OPARAM);
        return sum_col_counts_internal<Rcpp::IntegerVector>(mat.get(), out.get(), summation_set, start_index, end_index);

    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(counts);
        beachmat::output_param OPARAM(mat->get_class(), mat->get_package());
        auto out=beachmat::create_numeric_output(ngenes, nsummations, OPARAM);
        return sum_col_counts_internal<Rcpp::NumericVector>(mat.get(), out.get(), summation_set, start_index, end_index);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
} 

