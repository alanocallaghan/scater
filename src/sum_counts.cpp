#include "Rcpp.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/logical_matrix.h"

#include <stdexcept>

template <class M, class O>
Rcpp::RObject sum_row_counts_internal(Rcpp::RObject input, 
    const Rcpp::IntegerVector& genes, const Rcpp::IntegerVector& runs)
{
    auto mat=beachmat::create_matrix<M>(input);
    const size_t ncells=mat->get_ncol();
    typename M::vector holding(mat->get_nrow());

    const size_t nsummations=runs.size();
    O output(nsummations, ncells);

    for (size_t c=0; c<ncells; ++c) {
        mat->get_col(c, holding.begin());

        auto outcol=output.column(c);
        auto oIt=outcol.begin();
        
        auto it=holding.begin();
        auto gIt=genes.begin();

        for (auto s : runs) {
            while (s > 0) {
                (*oIt) += *(it + *gIt - 1); // get back to zero indexing.
                --s;
                ++gIt;
            }
            ++oIt;
        }
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::RObject sum_row_counts (Rcpp::RObject counts, Rcpp::IntegerVector genes, Rcpp::IntegerVector runs) {
    auto mattype=beachmat::find_sexp_type(counts);
    if (mattype==INTSXP) {
        return sum_row_counts_internal<beachmat::integer_matrix, Rcpp::IntegerMatrix>(counts, genes, runs);
    } else if (mattype==LGLSXP) {
        return sum_row_counts_internal<beachmat::logical_matrix, Rcpp::LogicalMatrix>(counts, genes, runs);
    } else if (mattype==REALSXP) {
        return sum_row_counts_internal<beachmat::numeric_matrix, Rcpp::NumericMatrix>(counts, genes, runs);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
} 
