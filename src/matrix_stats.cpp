#include "scater.h"

/* Function to calculate the variance across each row. */

template <typename T, class V, class M>
Rcpp::RObject row_vars_internal(M mat, Rcpp::RObject subset_row, Rcpp::RObject subset_col) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    Rcpp::IntegerVector rows=process_subset_vector(subset_row, mat, true);
    Rcpp::IntegerVector cols=process_subset_vector(subset_col, mat, false);
    Rcpp::NumericVector outvar(rows.size());
    
    if (cols.size() < 2) { 
        std::fill(outvar.begin(), outvar.end(), R_NaReal);
        return outvar;
    }
   
    V incoming(ncells);
    auto ovIt=outvar.begin();
    for (const auto& r : rows) {
        mat->get_row(r, incoming.begin());
        
        // Computing the mean.
        double curmean=0;
        for (const auto& c : cols) {
            curmean+=incoming[c];
        }
        curmean/=cols.size();

        // Computing the variance.
        double& curval=(*(ovIt++));
        for (const auto& c : cols) {
            const double tmp=incoming[c] - curmean;
            curval += tmp * tmp;
        }
        curval/=cols.size()-1;
    }

    return outvar;
}
 
SEXP row_vars(SEXP exprs, SEXP subset_row, SEXP subset_col) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return row_vars_internal<int, Rcpp::IntegerVector>(mat.get(), subset_row, subset_col);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return row_vars_internal<double, Rcpp::NumericVector>(mat.get(), subset_row, subset_col);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}
  
/* Function to calculate the variance across each column. */

template <typename T, class V, class M>
Rcpp::RObject col_vars_internal(M mat, Rcpp::RObject subset_row, Rcpp::RObject subset_col) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    Rcpp::IntegerVector rows=process_subset_vector(subset_row, mat, true);
    Rcpp::IntegerVector cols=process_subset_vector(subset_col, mat, false);
    Rcpp::NumericVector outvar(cols.size());
    
    if (rows.size() < 2) { 
        std::fill(outvar.begin(), outvar.end(), R_NaReal);
        return outvar;
    }

    V incoming(ngenes);
    auto ovIt=outvar.begin();
    for (const auto& c : cols) { 
        auto iIt=mat->get_const_col(c, incoming.begin());

        // Computing the mean.
        double curmean=0;
        for (const auto& r : rows) {
            curmean+=*(iIt+r);
        }
        curmean/=rows.size();

        // Computing the variance.
        double& curval=(*(ovIt++));
        for (const auto& r : rows) {
            const double tmp=*(iIt+r) - curmean;
            curval += tmp*tmp;
        }
        curval/=rows.size()-1;
    }

    return outvar;
}

SEXP col_vars(SEXP exprs, SEXP subset_row, SEXP subset_col) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return col_vars_internal<int, Rcpp::IntegerVector>(mat.get(), subset_row, subset_col);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return col_vars_internal<double, Rcpp::NumericVector>(mat.get(), subset_row, subset_col);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Function to calculate the sum across each row. */

template <typename T, class V, class M>
Rcpp::RObject row_sums_internal(M mat, Rcpp::RObject subset_row, Rcpp::RObject subset_col) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    Rcpp::IntegerVector rows=process_subset_vector(subset_row, mat, true);
    Rcpp::IntegerVector cols=process_subset_vector(subset_col, mat, false);
    V outsum(rows.size());
    if (!cols.size()) { 
        return outsum;
    }
   
    V incoming(ncells);
    auto osIt=outsum.begin();
    for (const auto& r : rows) {
        mat->get_row(r, incoming.begin());
        
        // Computing the mean.
        T& cursum=(*(osIt++));
        for (const auto& c : cols) {
            cursum+=incoming[c];
        }
    }

    return outsum;
}
 
SEXP row_sums(SEXP exprs, SEXP subset_row, SEXP subset_col) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return row_sums_internal<int, Rcpp::IntegerVector>(mat.get(), subset_row, subset_col);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return row_sums_internal<double, Rcpp::NumericVector>(mat.get(), subset_row, subset_col);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Function to calculate the sum across each column. */

template <typename T, class V, class M>
Rcpp::RObject col_sums_internal(M mat, Rcpp::RObject subset_row, Rcpp::RObject subset_col) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    Rcpp::IntegerVector rows=process_subset_vector(subset_row, mat, true);
    Rcpp::IntegerVector cols=process_subset_vector(subset_col, mat, false);
    V outsum(cols.size());
    if (!rows.size()) { 
        return outsum;
    }

    V incoming(ngenes);
    auto osIt=outsum.begin();
    for (const auto& c : cols) { 
        auto iIt=mat->get_const_col(c, incoming.begin());

        // Computing the mean.
        T& outsum=(*(osIt++));
        for (const auto& r : rows) {
            outsum+=*(iIt+r);
        }
    }

    return outsum;
}

SEXP col_sums(SEXP exprs, SEXP subset_row, SEXP subset_col) { 
    BEGIN_RCPP
    auto mattype=beachmat::find_sexp_type(exprs);
    if (mattype==INTSXP) {
        auto mat=beachmat::create_integer_matrix(exprs);
        return col_sums_internal<int, Rcpp::IntegerVector>(mat.get(), subset_row, subset_col);
    } else if (mattype==REALSXP) {
        auto mat=beachmat::create_numeric_matrix(exprs);
        return col_sums_internal<double, Rcpp::NumericVector>(mat.get(), subset_row, subset_col);
    } else {
        throw std::runtime_error("unacceptable matrix type");
    }
    END_RCPP
}

/* Function to count the incidences of particular value across each row. */

template <typename T, class V, class M>
Rcpp::RObject row_above_internal(M mat, Rcpp::RObject subset_row, Rcpp::RObject subset_col, Rcpp::RObject val) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    Rcpp::IntegerVector rows=process_subset_vector(subset_row, mat, true);
    Rcpp::IntegerVector cols=process_subset_vector(subset_col, mat, false);
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
   
    V incoming(ncells);
    auto ocIt=outcount.begin();
    for (const auto& r : rows) {
        mat->get_row(r, incoming.begin());
        
        int& curcount=(*(ocIt++));
        for (const auto& c : cols) {
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
Rcpp::RObject col_above_internal(M mat, Rcpp::RObject subset_row, Rcpp::RObject subset_col, Rcpp::RObject val) {
    const size_t& ncells=mat->get_ncol();
    const size_t& ngenes=mat->get_nrow();

    Rcpp::IntegerVector rows=process_subset_vector(subset_row, mat, true);
    Rcpp::IntegerVector cols=process_subset_vector(subset_col, mat, false);
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
 
    V incoming(ngenes);
    auto ocIt=outcount.begin();
    for (const auto& c : cols) { 
        auto iIt=mat->get_const_col(c, incoming.begin());

        int& curcount=(*(ocIt++));
        for (const auto& r : rows) {
            if (*(iIt+r) > objective) { 
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
