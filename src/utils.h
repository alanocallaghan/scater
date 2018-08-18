#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"

template<class MPTR>
Rcpp::IntegerVector process_subset_vector(Rcpp::RObject subset, MPTR mat, bool byrow) {
    if (subset.sexp_type()!=INTSXP) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }  
    Rcpp::IntegerVector sout(subset);
    const int& upper=(byrow ? mat->get_nrow() : mat->get_ncol());
    for (auto sIt=sout.begin(); sIt!=sout.end(); ++sIt) {
        if (*sIt < 0 || *sIt >= upper) { 
            throw std::runtime_error("subset indices out of range");
        }
    }
    return sout;
}

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

bool check_logical_scalar(Rcpp::RObject, const char*);

#endif

