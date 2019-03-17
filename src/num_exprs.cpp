#include "scater.h"

#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <utility>

/* Sorting rows and column subset indices for rapid matrix access. */

Rcpp::List reorganize_subset(Rcpp::IntegerVector sub1, Rcpp::IntegerVector sub2) {
    const size_t N=sub1.size();
    std::vector<std::pair<int, int> > paired(N);
    for (size_t i=0; i<N; ++i) {
        paired[i].first=sub1[i];
        paired[i].second=i;
    }
    std::sort(paired.begin(), paired.end()); // Not the fastest but duplicate-safe.

    Rcpp::IntegerVector sorted1(N), ordering(N);
    for (size_t i=0; i<N; ++i) {
        const auto& current=paired[i];
        sorted1[i]=current.first;
        ordering[i]=current.second;
    }

    Rcpp::IntegerVector sorted2 = Rcpp::clone(sub2).sort();
    int first=sorted2[0], last=sorted2[sorted2.size()-1]+1; // sub2 should be non-empty by this point.
    for (auto& x : sorted2) { 
        x-=first;
    }

    return Rcpp::List::create(sorted1, ordering, sorted2, 
            Rcpp::IntegerVector::create(first, last));
}

/* Function to count the incidences of particular value across each row. */

template <typename T, class V, class M>
Rcpp::RObject row_above_internal(M mat, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols, Rcpp::RObject val) {
    Rcpp::IntegerVector outcount(rows.size());
    if (!cols.size()) { 
        return outcount;
    }

    // Coercing the target to the same type as the matrix.
    V target(val);
    if (target.size()!=1) { 
        throw std::runtime_error("value to find must be a scalar");
    }
    const T& objective=target[0];

    Rcpp::List reorganized=reorganize_subset(rows, cols);
    Rcpp::IntegerVector rowcopy(reorganized[0]);
    Rcpp::IntegerVector roworder(reorganized[1]);
    Rcpp::IntegerVector colcopy(reorganized[2]);
    Rcpp::IntegerVector offsets(reorganized[3]);
    int first=offsets[0], last=offsets[1]; 

    V incoming(mat->get_ncol());
    auto roIt=roworder.begin();
    for (const auto& r : rowcopy) {
        mat->get_row(r, incoming.begin(), first, last);
        
        int& curcount=outcount[*(roIt++)];
        for (const auto& c : colcopy) {
            if (incoming[c] > objective) { 
                ++curcount;
            }
        }
    }

    return outcount;
}
 
SEXP row_above(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP value) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return row_above_internal<int, Rcpp::IntegerVector>(mat.get(), subset_row, subset_col, value);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return row_above_internal<double, Rcpp::NumericVector>(mat.get(), subset_row, subset_col, value);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Function to count the incidences of particular value across each column. */

template <typename T, class V, class M>
Rcpp::RObject col_above_internal(M mat, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols, Rcpp::RObject val) {
    Rcpp::IntegerVector outcount(cols.size());
    if (!rows.size()) { 
        return outcount;
    }

    // Coercing the target to the same type as the matrix.
    V target(val);
    if (target.size()!=1) { 
        throw std::runtime_error("value to find must be a scalar");
    }
    const T& objective=target[0];

    Rcpp::List reorganized=reorganize_subset(cols, rows);
    Rcpp::IntegerVector colcopy(reorganized[0]);
    Rcpp::IntegerVector colorder(reorganized[1]);
    Rcpp::IntegerVector rowcopy(reorganized[2]);
    Rcpp::IntegerVector offsets(reorganized[3]);
    int first=offsets[0], last=offsets[1]; 
 
    V incoming(mat->get_nrow());
    auto coIt=colorder.begin();
    for (const auto& c : colcopy) { 
        mat->get_col(c, incoming.begin(), first, last);

        int& curcount=outcount[*(coIt++)];
        for (const auto& r : rowcopy) {
            if (incoming[r] > objective) { 
                ++curcount;
            }
        }
    }

    return outcount;
}

SEXP col_above(SEXP exprs, SEXP subset_row, SEXP subset_col, SEXP value) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return col_above_internal<int, Rcpp::IntegerVector>(mat.get(), subset_row, subset_col, value);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return col_above_internal<double, Rcpp::NumericVector>(mat.get(), subset_row, subset_col, value);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
