#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"

Rcpp::IntegerVector process_subset_vector(Rcpp::RObject, int, bool=true);

int check_integer_scalar(Rcpp::RObject, const char*);

double check_numeric_scalar(Rcpp::RObject, const char*);

bool check_logical_scalar(Rcpp::RObject, const char*);

#endif

