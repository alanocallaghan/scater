#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"
#include "beachmat/utils/const_column.h"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <utility>

/* Function to count the incidences of particular value across each row. */

template <class M>
Rcpp::RObject row_above_internal(Rcpp::RObject input, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols, Rcpp::RObject val) {
    auto mat=beachmat::create_matrix<M>(input);
    Rcpp::IntegerVector outcount(rows.size());
    if (!cols.size()) { 
        return outcount;
    }

    // Coercing the target to the same type as the matrix.
    typename M::vector target(val);
    if (target.size()!=1) { 
        throw std::runtime_error("value to find must be a scalar");
    }
    typename M::type objective=target[0];
 
    int first=0, last=0;
    if (rows.size()) {
        first=*std::min_element(rows.begin(), rows.end());
        last=*std::max_element(rows.begin(), rows.end())+1;
    }

    // Difficult to handle subsetting with const column getters,
    // so we disable it here.
    beachmat::const_column<M> col_holder(mat.get(), false);

    for (const auto c : cols) {
        col_holder.fill(c, first, last);
        auto it=col_holder.get_values();
        auto oIt=outcount.begin();
        for (const auto r : rows) {
            if (*(it + r - first) > objective) { ++(*oIt); }
            ++oIt;
        }
    }

    return outcount;
}
 
SEXP row_above(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP value) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        return row_above_internal<beachmat::integer_matrix>(exprs, subset_row, subset_col, value);
    } else if (mattype==REALSXP) {
        return row_above_internal<beachmat::numeric_matrix>(exprs, subset_row, subset_col, value);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Function to count the incidences of particular value across each column. */

template <class M>
Rcpp::RObject col_above_internal(Rcpp::RObject input, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols, Rcpp::RObject val) {
    auto mat=beachmat::create_matrix<M>(input);
    Rcpp::IntegerVector outcount(cols.size());
    if (!rows.size()) { 
        return outcount;
    }

    // Coercing the target to the same type as the matrix.
    typename M::vector target(val);
    if (target.size()!=1) { 
        throw std::runtime_error("value to find must be a scalar");
    }
    typename M::type objective=target[0];

    int first=0, last=0;
    if (rows.size()) {
        first=*std::min_element(rows.begin(), rows.end());
        last=*std::max_element(rows.begin(), rows.end())+1;
    }

    beachmat::const_column<M> col_holder(mat.get(), false);
    auto oIt=outcount.begin();
    for (const auto c : cols) { 
        col_holder.fill(c, first, last);
        auto it=col_holder.get_values();
        for (const auto r : rows) {
            if (*(it + r - first) > objective) { ++(*oIt); }
        }
        ++oIt;
    }

    return outcount;
}

SEXP col_above(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP value) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        return col_above_internal<beachmat::integer_matrix>(exprs, subset_row, subset_col, value);
    } else if (mattype==REALSXP) {
        return col_above_internal<beachmat::numeric_matrix>(exprs, subset_row, subset_col, value);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
