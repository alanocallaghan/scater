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
#include <algorithm>
#include <stdexcept>

SEXP Make(R_altrep_class_t* class_t, SEXP mat, SEXP dim, SEXP idx) {
    MARK_NOT_MUTABLE(mat);
    MARK_NOT_MUTABLE(dim);
    MARK_NOT_MUTABLE(idx);

    SEXP data1=PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(data1, 0, mat);
    SET_VECTOR_ELT(data1, 1, dim);
    SET_VECTOR_ELT(data1, 2, idx);

    SEXP out=R_new_altrep(*class_t, data1, R_NilValue);
    UNPROTECT(1);
    return out;
}

template<typename GRABBER>
struct lazy_vector_methods {
    static SEXP Materialize(SEXP vec) {
        SEXP data2=R_altrep_data2(vec);
        if (data2==R_NilValue) {
            SEXP data1=R_altrep_data1(vec);
            GRABBER grab(data1);

            int veclen=grab.size();
            data2=PROTECT(Rf_allocVector(GRABBER::sexp_type, veclen));
            grab.fill(static_cast<typename GRABBER::value_type*>(DATAPTR(data2)));

            R_set_altrep_data2(vec, data2);
            UNPROTECT(1);
        }
        return data2;
    }

    // ALTREP methods -------------------

    static R_xlen_t Length(SEXP vec) {
        SEXP data2=R_altrep_data2(vec);
        if (data2==R_NilValue) {
            SEXP data1=R_altrep_data1(vec);
            GRABBER grab(data1);
            return grab.size();
        } else {
            return LENGTH(data2);
        }
    }

    static Rboolean Inspect(SEXP x, int pre, int deep, int pvec, void (*inspect_subtree)(SEXP, int, int, int)){
        auto len=LENGTH(x);
        const char* dimname=(GRABBER::getcol ? "row" : "column");
        const char* typname=(GRABBER::sexp_type==REALSXP ? "double" : "integer");

        SEXP data2=R_altrep_data2(x);
        const char* mode=(data2==R_NilValue ? "lazy" : "materialized");
        Rprintf("%s %s (len=%d, type=%s)\n", mode, dimname, len, typname);

        return TRUE;
    }

    // ALTVEC methods ------------------

    static const void* Dataptr_or_null(SEXP vec){
        SEXP data2 = R_altrep_data2(vec);
        if (data2 == R_NilValue) {
            return nullptr;
        } else {
            return STDVEC_DATAPTR(data2);
        }
    }

    static void* Dataptr(SEXP vec, Rboolean writeable){
        return STDVEC_DATAPTR(Materialize(vec));
    }

    // ALTREAL/INTEGER methods -----------------

    static typename GRABBER::value_type value_Elt(SEXP vec, R_xlen_t i){
        SEXP data2=R_altrep_data2(vec);
        if (data2==R_NilValue) {
            SEXP data1=R_altrep_data1(vec);
            GRABBER grab(data1);
            return grab(i);
        } else {
            return static_cast<const typename GRABBER::value_type*>(STDVEC_DATAPTR(data2))[i];
        }
    }

    static R_xlen_t Get_region(SEXP vec, R_xlen_t start, R_xlen_t size, typename GRABBER::value_type* out){
        out = static_cast<typename GRABBER::value_type*>(Dataptr(vec, TRUE));
        R_xlen_t len = Length(vec) - start;
        return len > size ? len : size;
    }

    static void partial_Init(R_altrep_class_t* class_t) {
        R_set_altrep_Length_method(*class_t, Length);
        R_set_altrep_Inspect_method(*class_t, Inspect);
        R_set_altvec_Dataptr_method(*class_t, Dataptr);
        R_set_altvec_Dataptr_or_null_method(*class_t, Dataptr_or_null);
    }
};

template<typename METHODS>
void remaining_int_Init(R_altrep_class_t* class_t) {
    R_set_altinteger_Elt_method(*class_t, METHODS::value_Elt);
    R_set_altinteger_Get_region_method(*class_t, METHODS::Get_region);
}

template<typename METHODS>
void remaining_dbl_Init(R_altrep_class_t* class_t) {
    R_set_altreal_Elt_method(*class_t, METHODS::value_Elt);
    R_set_altreal_Get_region_method(*class_t, METHODS::Get_region);
}

/***************************************************/

template<bool GETCOL>
struct base_grabber {
    static const bool getcol=GETCOL;

    base_grabber(SEXP data1) {
        raw_mat=VECTOR_ELT(data1, 0);
        dim=VECTOR_ELT(data1, 1);
        idx=Rcpp::as<int>(VECTOR_ELT(data1, 2));
        return;
    }

    int size() const {
        // Flipped around, because if we want the column,
        // then obviously we need the number of rows.
        if (GETCOL) {
            return dim[0];
        } else {
            return dim[1];
        }
    }
protected:
    Rcpp::RObject raw_mat;
    Rcpp::IntegerVector dim;
    int idx;
};


template<bool GETCOL, typename M>
struct ordinary_grabber : public base_grabber<GETCOL> {
    typedef typename M::stored_type value_type;
    static const unsigned int sexp_type=M::r_type::value;

    ordinary_grabber(SEXP data1) : base_grabber<GETCOL>(data1), mat(this->raw_mat) {}

    value_type operator()(int i) {
        if (GETCOL) {
            return mat(i, this->idx);
        } else {
            return mat(this->idx, i);
        }
    }

    void fill(value_type* out) {
        if (GETCOL) {
            auto C=mat.column(this->idx);
            std::copy(C.begin(), C.end(), out);
        } else {
            auto R=mat.row(this->idx);
            std::copy(R.begin(), R.end(), out);
        }
        return;
    }
protected:
    M mat;
};

template<bool GETCOL, typename V>
struct sparse_grabber : base_grabber<GETCOL> {
    typedef typename V::stored_type value_type;
    static const unsigned int sexp_type=V::r_type::value;

    sparse_grabber(SEXP data1) : base_grabber<GETCOL>(data1) {
        Rcpp::S4 tmp(this->raw_mat);
        x_=tmp.slot("x");
        i_=tmp.slot("i");
        p_=tmp.slot("p");
        return;
    }

    value_type operator()(int i) {
        if (GETCOL) {
            return search(i, this->idx);
        } else {
            return search(this->idx, i);
        }
    }

    void fill(value_type* out) {
        if (GETCOL) {
            std::fill(out, out+this->dim[0], 0);
            auto xIt=x_.begin();
            auto iIt=i_.begin();
            for (int x=p_[this->idx]; x!=p_[this->idx+1]; ++x) {
                out[*(iIt+x)]=*(xIt+x);
            }
        } else {
            for (int j=0; j<this->dim[1]; ++j) {
                out[j]=search(this->idx, j);
            }
        }
        return;
    }
protected:
    V x_;
    Rcpp::IntegerVector i_;
    Rcpp::IntegerVector p_;

    double search(int i, int j) {
        int start=p_[j], end=p_[j+1];
        auto last=i_.begin() + end;
        auto out=std::lower_bound(i_.begin() + start, last, i);
        if (out==last || *out != i) {
            return 0;
        } else {
            return x_[out - i_.begin()];
        }
    }
};

/***************************************************/

typedef ordinary_grabber<true, Rcpp::IntegerMatrix> ordinary_int_col_grabber;

typedef ordinary_grabber<false, Rcpp::IntegerMatrix> ordinary_int_row_grabber;

typedef ordinary_grabber<true, Rcpp::NumericMatrix> ordinary_dbl_col_grabber;

typedef ordinary_grabber<false, Rcpp::NumericMatrix> ordinary_dbl_row_grabber;

typedef sparse_grabber<true, Rcpp::NumericVector> sparse_dbl_col_grabber;

typedef sparse_grabber<false, Rcpp::NumericVector> sparse_dbl_row_grabber;

/***************************************************/

typedef lazy_vector_methods<ordinary_int_col_grabber> lazy_ordinary_int_col_methods;

typedef lazy_vector_methods<ordinary_int_row_grabber> lazy_ordinary_int_row_methods;

typedef lazy_vector_methods<ordinary_dbl_col_grabber> lazy_ordinary_dbl_col_methods;

typedef lazy_vector_methods<ordinary_dbl_row_grabber> lazy_ordinary_dbl_row_methods;

typedef lazy_vector_methods<sparse_dbl_col_grabber> lazy_sparse_dbl_col_methods;

typedef lazy_vector_methods<sparse_dbl_row_grabber> lazy_sparse_dbl_row_methods;

/***************************************************/

R_altrep_class_t lazy_ordinary_int_col_t;
R_altrep_class_t lazy_ordinary_int_row_t;
R_altrep_class_t lazy_ordinary_dbl_col_t;
R_altrep_class_t lazy_ordinary_dbl_row_t;
R_altrep_class_t lazy_sparse_dbl_col_t;
R_altrep_class_t lazy_sparse_dbl_row_t;

// [[Rcpp::init]]
void init_lazy_vector(DllInfo* dll){
    lazy_ordinary_int_col_t = R_make_altinteger_class("lazy_ordinary_int_col", "scater", dll);
    lazy_ordinary_int_col_methods::partial_Init(&lazy_ordinary_int_col_t);
    remaining_int_Init<lazy_ordinary_int_col_methods>(&lazy_ordinary_int_col_t);

    lazy_ordinary_int_row_t = R_make_altinteger_class("lazy_ordinary_int_row", "scater", dll);
    lazy_ordinary_int_row_methods::partial_Init(&lazy_ordinary_int_row_t);
    remaining_int_Init<lazy_ordinary_int_row_methods>(&lazy_ordinary_int_row_t);

    lazy_ordinary_dbl_col_t = R_make_altreal_class("lazy_ordinary_dbl_col", "scater", dll);
    lazy_ordinary_dbl_col_methods::partial_Init(&lazy_ordinary_dbl_col_t);
    remaining_dbl_Init<lazy_ordinary_dbl_col_methods>(&lazy_ordinary_dbl_col_t);

    lazy_ordinary_dbl_row_t = R_make_altreal_class("lazy_ordinary_dbl_row", "scater", dll);
    lazy_ordinary_dbl_row_methods::partial_Init(&lazy_ordinary_dbl_row_t);
    remaining_dbl_Init<lazy_ordinary_dbl_row_methods>(&lazy_ordinary_dbl_row_t);

    lazy_sparse_dbl_col_t = R_make_altreal_class("lazy_sparse_dbl_col", "scater", dll);
    lazy_sparse_dbl_col_methods::partial_Init(&lazy_sparse_dbl_col_t);
    remaining_dbl_Init<lazy_sparse_dbl_col_methods>(&lazy_sparse_dbl_col_t);

    lazy_sparse_dbl_row_t = R_make_altreal_class("lazy_sparse_dbl_row", "scater", dll);
    lazy_sparse_dbl_row_methods::partial_Init(&lazy_sparse_dbl_row_t);
    remaining_dbl_Init<lazy_sparse_dbl_row_methods>(&lazy_sparse_dbl_row_t);
}

// [[Rcpp::export(rng=false)]]
SEXP create_lazy_vector(SEXP mat, SEXP dim, SEXP idx, bool getcol, int matclass, int type) {
    /*
     * matclass = 0 is an ordinary matrix.
     * matclass = 1 is a dgCMatrix.
     *
     * type = 0 is integer.
     * type = 1 is double.
     */

    if (matclass==0) {
        if (type==0) {
            if (getcol) {
                return Make(&lazy_ordinary_int_col_t, mat, dim, idx);
            } else {
                return Make(&lazy_ordinary_int_row_t, mat, dim, idx);
            }
        } else if (type==1) {
            if (getcol) {
                return Make(&lazy_ordinary_dbl_col_t, mat, dim, idx);
            } else {
                return Make(&lazy_ordinary_dbl_row_t, mat, dim, idx);
            }
        }
    } else if (matclass==1) {
        if (type==1) {
            if (getcol) {
                return Make(&lazy_sparse_dbl_col_t, mat, dim, idx);
            } else {
                return Make(&lazy_sparse_dbl_row_t, mat, dim, idx);
            }
        }
    }
    throw std::runtime_error("lazy vectors not supported for this assay");
}
