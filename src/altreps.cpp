#include "Rcpp.h"
#include "R.h"
#include "Rdefines.h"
#include <R_ext/Rdynload.h>

#define class klass

extern "C" {
#include "R_ext/Altrep.h"
}

#undef klass

#include <type_traits>

SEXP Make(R_altrep_class_t* class_t, SEXP mat, SEXP idx) {
    MARK_NOT_MUTABLE(mat);
    MARK_NOT_MUTABLE(idx);

    SEXP data1=PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(data1, 0, mat);
    SET_VECTOR_ELT(data1, 1, idx);

    SEXP out=R_new_altrep(*class_t, data1, R_NilValue);
    UNPROTECT(1);
    return out;
}

template<typename T, bool GET_COL>
SEXP Materialize(SEXP vec) {
    SEXP data2=R_altrep_data2(vec);
    if (data2==R_NilValue) {
        SEXP data1=R_altrep_data1(vec);
        SEXP mat=VECTOR_ELT(data1, 0);
        SEXP idx=VECTOR_ELT(data1, 1);

        Rcpp::Environment env=Rcpp::Environment::namespace_env("scater");
        Rcpp::RObject realized;
        if (GET_COL) {
            Rcpp::Function FUN(env["getColumn"]);
            realized=FUN(mat, idx);
        } else {
            Rcpp::Function FUN(env["getRow"]);
            realized=FUN(mat, idx);
        }

        R_set_altrep_data2(vec, SEXP(realized));
        data2=R_altrep_data2(vec);
    }
    return data2;
}

// ALTREP methods -------------------

template<bool GET_COL>
static R_xlen_t Length(SEXP vec) {
    SEXP data2=R_altrep_data2(vec);
    if (data2==R_NilValue) {
        SEXP data1=R_altrep_data1(vec);
        SEXP mat=VECTOR_ELT(data1, 0);
        Rcpp::Environment env=Rcpp::Environment::namespace_env("base");
        Rcpp::RObject output;

        const char* funname=(GET_COL ? "nrow" : "ncol");
        Rcpp::Function FUN(env[funname]);
        output=FUN(mat);
        return Rcpp::as<int>(output);
    } else {
        return LENGTH(data2);
    }
}

template<typename T, bool GET_COL>
Rboolean Inspect(SEXP x, int pre, int deep, int pvec, void (*inspect_subtree)(SEXP, int, int, int)){
    auto len=LENGTH(x);
    const char* dimname=(GET_COL ? "row" : "column");
    const char* typname=(std::is_same<T, double>::value ? "double" : "integer");

    SEXP data2=R_altrep_data2(x);
    const char* mode=(data2==R_NilValue ? "lazy" : "materialized");
    Rprintf("%s %s (len=%d, type=%s)\n", mode, dimname, len, typname);

    return TRUE;
}

// ALTVEC methods ------------------

const void* Dataptr_or_null(SEXP vec){
    SEXP data2 = R_altrep_data2(vec);
    if (data2 == R_NilValue) {
        return nullptr;
    } else {
        return STDVEC_DATAPTR(data2);
    }
}

template<typename T, bool GET_COL>
void* Dataptr(SEXP vec, Rboolean writeable){
    return STDVEC_DATAPTR(Materialize<T, GET_COL>(vec));
}

// ALTREAL/INTEGER methods -----------------

template<typename T, bool GET_COL>
static T value_Elt(SEXP vec, R_xlen_t i){
    SEXP data2=R_altrep_data2(vec);
    if (data2==R_NilValue) {
        data2=Materialize<T, GET_COL>(vec);
    } 
    return static_cast<T*>(STDVEC_DATAPTR(data2))[i];
}

template<typename T, bool GET_COL>
static R_xlen_t Get_region(SEXP vec, R_xlen_t start, R_xlen_t size, T* out){
    T* data = static_cast<T*>(Dataptr<T, GET_COL>(vec, TRUE));
    out = data + start;
    R_xlen_t len = Length<GET_COL>(vec) - start;
    return len > size ? len : size;
}

/***************************************************/

R_altrep_class_t int_col_t;
R_altrep_class_t int_row_t;
R_altrep_class_t dbl_col_t;
R_altrep_class_t dbl_row_t;

// [[Rcpp::init]]
void init_lazy_vector(DllInfo* dll){
    int_col_t = R_make_altinteger_class("lazy_integer_column", "scater", dll);
    R_set_altrep_Length_method(int_col_t, Length<true>);
    R_set_altrep_Inspect_method(int_col_t, Inspect<int, true>);
    R_set_altvec_Dataptr_method(int_col_t, Dataptr<int, true>);
    R_set_altvec_Dataptr_or_null_method(int_col_t, Dataptr_or_null);
    R_set_altinteger_Elt_method(int_col_t, value_Elt<int, true>);
    R_set_altinteger_Get_region_method(int_col_t, Get_region<int, true>);

    int_row_t = R_make_altinteger_class("lazy_integer_row", "scater", dll);
    R_set_altrep_Length_method(int_row_t, Length<false>);
    R_set_altrep_Inspect_method(int_row_t, Inspect<int, false>);
    R_set_altvec_Dataptr_method(int_row_t, Dataptr<int, false>);
    R_set_altvec_Dataptr_or_null_method(int_row_t, Dataptr_or_null);
    R_set_altinteger_Elt_method(int_row_t, value_Elt<int, false>);
    R_set_altinteger_Get_region_method(int_row_t, Get_region<int, false>);

    dbl_col_t = R_make_altreal_class("lazy_double_column", "scater", dll);
    R_set_altrep_Length_method(dbl_col_t, Length<true>);
    R_set_altrep_Inspect_method(dbl_col_t, Inspect<double, true>);
    R_set_altvec_Dataptr_method(dbl_col_t, Dataptr<double, true>);
    R_set_altvec_Dataptr_or_null_method(dbl_col_t, Dataptr_or_null);
    R_set_altreal_Elt_method(dbl_col_t, value_Elt<double, true>);
    R_set_altreal_Get_region_method(dbl_col_t, Get_region<double, true>);

    dbl_row_t = R_make_altreal_class("lazy_double_row", "scater", dll);
    R_set_altrep_Length_method(dbl_row_t, Length<false>);
    R_set_altrep_Inspect_method(dbl_row_t, Inspect<double, false>);
    R_set_altvec_Dataptr_method(dbl_row_t, Dataptr<double, false>);
    R_set_altvec_Dataptr_or_null_method(dbl_row_t, Dataptr_or_null);
    R_set_altreal_Elt_method(dbl_row_t, value_Elt<double, false>);
    R_set_altreal_Get_region_method(dbl_row_t, Get_region<double, false>);
}

// [[Rcpp::export(rng=false)]]
SEXP lazy_integer_column(SEXP mat, SEXP idx) {
    return Make(&int_col_t, mat, idx);
}

// [[Rcpp::export(rng=false)]]
SEXP lazy_integer_row(SEXP mat, SEXP idx) {
    return Make(&int_row_t, mat, idx);
}

// [[Rcpp::export(rng=false)]]
SEXP lazy_double_column(SEXP mat, SEXP idx) {
    return Make(&dbl_col_t, mat, idx);
}

// [[Rcpp::export(rng=false)]]
SEXP lazy_double_row(SEXP mat, SEXP idx) {
    return Make(&dbl_row_t, mat, idx);
}
