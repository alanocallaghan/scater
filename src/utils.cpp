#include "utils.h"

// Checking subsetting vector.

Rcpp::IntegerVector process_subset_vector(Rcpp::RObject subset, int upper, bool zero_indexed) {
    if (subset.sexp_type()!=INTSXP) { 
        throw std::runtime_error("subset vector must be an integer vector");
    }  
    Rcpp::IntegerVector sout(subset);

    if (!zero_indexed) {
        sout=Rcpp::clone(sout);
        for (auto& s : sout) { --s; }
    }

    for (auto sIt=sout.begin(); sIt!=sout.end(); ++sIt) {
        if (*sIt < 0 || *sIt >= upper) { 
            throw std::runtime_error("subset indices out of range");
        }
    }
    return sout;
}

// Checking for scalar inputs.

template<typename T, class V>
T check_scalar(Rcpp::RObject incoming, const char* arg, const char* val) {
    V vec(incoming);
    if (vec.size()!=1) {
        std::stringstream err;
        err << arg << " should be " << val;
        throw std::runtime_error(err.str());
    }
    return vec[0];
}

int check_integer_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<int, Rcpp::IntegerVector>(incoming, arg, "an integer scalar");
}

double check_numeric_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<double, Rcpp::NumericVector>(incoming, arg, "a numeric scalar");
}

bool check_logical_scalar(Rcpp::RObject incoming, const char* arg) {
    return check_scalar<bool, Rcpp::LogicalVector>(incoming, arg, "a logical scalar");
}
