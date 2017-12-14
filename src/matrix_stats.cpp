#include "scater.h"

/* Sorting rows and column subset indices for rapid matrix access. */

Rcpp::List reorganize_subset(Rcpp::IntegerVector sub1, Rcpp::IntegerVector sub2) {
    Rcpp::IntegerVector sorted1 = Rcpp::clone(sub1).sort();
    Rcpp::IntegerVector ordering = Rcpp::match(sorted1, sub1);  
    for (auto& o : ordering) { 
        --o;
    }

    Rcpp::IntegerVector sorted2 = Rcpp::clone(sub2).sort();
    int first=sorted2[0], last=sorted2[sorted2.size()-1]+1; // sub2 should be non-empty by this point.
    for (auto& x : sorted2) { 
        x-=first;
    }

    return Rcpp::List::create(sorted1, ordering, sorted2, 
            Rcpp::IntegerVector::create(first, last));
}

/* Function to calculate the variance across each row. */

template <typename T, class V, class M>
Rcpp::RObject row_vars_internal(M mat, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols) {
    if (cols.size() < 2) { 
        return Rcpp::NumericVector(rows.size(), R_NaReal);
    }

    Rcpp::List reorganized=reorganize_subset(rows, cols);
    Rcpp::IntegerVector rowcopy(reorganized[0]);
    Rcpp::IntegerVector roworder(reorganized[1]);
    Rcpp::IntegerVector colcopy(reorganized[2]);
    Rcpp::IntegerVector offsets(reorganized[3]);
    int first=offsets[0], last=offsets[1]; 

    Rcpp::NumericVector outvar(rows.size());  
    V incoming(mat->get_ncol());
    auto roIt=roworder.begin();

    for (const auto& r : rowcopy) {
        mat->get_row(r, incoming.begin(), first, last);
        
        // Computing the mean.
        double curmean=0;
        for (const auto& c : colcopy) {
            curmean+=incoming[c];
        }
        curmean/=cols.size();

        // Computing the variance.
        double& curval=outvar[*(roIt++)];
        for (const auto& c : colcopy) {
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
Rcpp::RObject col_vars_internal(M mat, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols) { 
    if (rows.size() < 2) { 
        return Rcpp::NumericVector(cols.size(), R_NaReal);
    }
   
    Rcpp::List reorganized=reorganize_subset(cols, rows);
    Rcpp::IntegerVector colcopy(reorganized[0]);
    Rcpp::IntegerVector colorder(reorganized[1]);
    Rcpp::IntegerVector rowcopy(reorganized[2]);
    Rcpp::IntegerVector offsets(reorganized[3]);
    int first=offsets[0], last=offsets[1]; 
 
    Rcpp::NumericVector outvar(cols.size()); 
    V incoming(mat->get_nrow());
    auto coIt=colorder.begin();
    for (const auto& c : colcopy) { 
        auto iIt=mat->get_const_col(c, incoming.begin(), first, last);

        // Computing the mean.
        double curmean=0;
        for (const auto& r : rowcopy) {
            curmean+=*(iIt+r);
        }
        curmean/=rows.size();

        // Computing the variance.
        double& curval=outvar[*(coIt++)];
        for (const auto& r : rowcopy) {
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
Rcpp::RObject row_sums_internal(M mat, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols) {
    V outsum(rows.size());
    if (!cols.size()) { 
        return outsum;
    }

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
        
        // Computing the mean.
        T& cursum=outsum[*(roIt++)];
        for (const auto& c : colcopy) {
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
Rcpp::RObject col_sums_internal(M mat, Rcpp::IntegerVector rows, Rcpp::IntegerVector cols) { 
    V outsum(cols.size());
    if (!rows.size()) { 
        return outsum;
    }
    
    Rcpp::List reorganized=reorganize_subset(cols, rows);
    Rcpp::IntegerVector colcopy(reorganized[0]);
    Rcpp::IntegerVector colorder(reorganized[1]);
    Rcpp::IntegerVector rowcopy(reorganized[2]);
    Rcpp::IntegerVector offsets(reorganized[3]);
    int first=offsets[0], last=offsets[1]; 
 
    V incoming(mat->get_nrow());
    auto coIt=colorder.begin();
    for (const auto& c : colcopy) { 
        auto iIt=mat->get_const_col(c, incoming.begin(), first, last);

        // Computing the mean.
        T& cursum=outsum[*(coIt++)];
        for (const auto& r : rowcopy) {
            cursum+=*(iIt+r);
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
        auto iIt=mat->get_const_col(c, incoming.begin(), first, last);

        int& curcount=outcount[*(coIt++)];
        for (const auto& r : rowcopy) {
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
