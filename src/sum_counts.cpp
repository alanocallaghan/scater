#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"
#include "utils.h"

#include <stdexcept>
#include <vector>

/***************************
 * Summing across features *
 ***************************/

template <class M, class O>
Rcpp::RObject sum_row_counts_internal(Rcpp::RObject input, 
    const Rcpp::IntegerVector& genes, const Rcpp::IntegerVector& runs,
    size_t start_index, size_t end_index) 
{
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ncells=mat->get_ncol();

    const size_t nsummations=runs.size();
    typename M::vector holder_out(nsummations);
    if (end_index > ncells) {
        throw std::runtime_error("end index out of range");
    }

    beachmat::output_param OPARAM(mat.get());
    const size_t njobs=end_index - start_index;
    auto out=beachmat::create_output<O>(nsummations, njobs, OPARAM);

    // Specializing for dense arrays to avoid copies.
    // Again, quite difficult to support sparsity due to the
    // need for arbitrary indexing of rows.
    beachmat::const_column<M> col_holder(mat.get(), false);

    for (size_t c=start_index; c<end_index; ++c) {
        col_holder.fill(c);

        auto it=col_holder.get_values();
        auto hIt=holder_out.begin();
        auto gIt=genes.begin();

        for (auto s : runs) {
            while (s > 0) {
                (*hIt) += *(it + *gIt - 1); // get back to zero indexing.
                --s;
                ++gIt;
            }
            ++hIt;
        }

        out->set_col(c - start_index, holder_out.begin());
        std::fill(holder_out.begin(), holder_out.end(), 0);
    }
    return out->yield();
}

SEXP sum_row_counts (SEXP counts, SEXP genes, SEXP runs, SEXP job_start, SEXP job_end) {
    BEGIN_RCPP

    // Determining the type and amount of work to do.
    size_t start_index=check_integer_scalar(job_start, "start index");
    size_t end_index=check_integer_scalar(job_end, "end index");
    if (end_index < start_index) { 
        throw std::runtime_error("start index is less than end index"); 
    }

    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        return sum_row_counts_internal<beachmat::integer_matrix,
           beachmat::integer_output>(counts, genes, runs, start_index, end_index);

    } else if (mattype==REALSXP) {
        return sum_row_counts_internal<beachmat::numeric_matrix,
           beachmat::numeric_output>(counts, genes, runs, start_index, end_index);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
} 

/************************
 * Summing across cells *
 ************************/

template <class M, class O>
Rcpp::RObject sum_col_counts_internal(Rcpp::RObject input, const Rcpp::List& summable_set, 
    size_t start_index, size_t end_index)
{
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ncells=mat->get_ncol();
    const size_t ngenes=mat->get_nrow();

    typename M::vector holder_out(ngenes);
    if (end_index > ngenes) {
        throw std::runtime_error("end index out of range");
    }
    beachmat::const_column<M> col_holder(mat.get());

    beachmat::output_param OPARAM(mat.get());
    const size_t njobs=end_index - start_index;
    const size_t nsummations=summable_set.size();
    auto out=beachmat::create_output<O>(njobs, nsummations, OPARAM);

    /* We do the rather unintuitive approach of running across each set
     * and summing columns, rather than running across columns and adding
     * them to the different sets as we go along. This approach has the
     * benefit of avoiding the need to hold intermediate sums, which may
     * be costly to re-extract depending on the backend used in output_param.
     */
    for (size_t i_summed=0; i_summed!=summable_set.size(); ++i_summed) {
        Rcpp::IntegerVector curset=summable_set[i_summed];
        std::fill(holder_out.begin()+start_index, holder_out.begin()+end_index, 0);

        for (auto cell : curset) {
            col_holder.fill(cell, start_index, end_index);
            auto val=col_holder.get_values();

            // Specialized for sparse matrices.
            if (col_holder.is_sparse()) {
                auto n=col_holder.get_n();
                auto idx=col_holder.get_indices();
                for (size_t i=0; i<n; ++i, ++idx, ++val) {
                    holder_out[*idx]+=*val;
                }
            } else {
                for (size_t i=start_index; i<end_index; ++i, ++val) {
                    holder_out[i] += *val;
                }
            }
        } 

        out->set_col(i_summed, holder_out.begin()+start_index);
    }
    return out->yield();
}

SEXP sum_col_counts (SEXP counts, SEXP sumset, SEXP job_start, SEXP job_end) {
    BEGIN_RCPP

    // Determining the type and amount of work to do.
    size_t start_index=check_integer_scalar(job_start, "start index");
    size_t end_index=check_integer_scalar(job_end, "end index");
    if (end_index < start_index) { 
        throw std::runtime_error("start index is less than end index"); 
    }

    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        return sum_col_counts_internal<beachmat::integer_matrix,
           beachmat::integer_output>(counts, sumset, start_index, end_index);

    } else if (mattype==REALSXP) {
        return sum_col_counts_internal<beachmat::numeric_matrix,
           beachmat::numeric_output>(counts, sumset, start_index, end_index);

    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
} 
