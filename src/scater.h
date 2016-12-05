#ifndef SCATER_H
#define SCATER_H

#include "R.h"
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"
#include "Rmath.h"

#include <algorithm>
#include <deque>
#include <stdexcept>
#include <functional>
#include <cmath>

// Facilitate matrix input.

struct matrix_info {
    matrix_info(int, int, bool);
    const size_t nrow, ncol;
    const bool is_integer;
    const int* iptr;
    const double* dptr;
};

matrix_info check_matrix(SEXP matrix);

// Special function to check for NA'ness.

bool isNA(int);
bool isNA(double);

// Functions to parse subsetting.

typedef std::pair<const int, const int*> subset_info;
subset_info process_subset_vector(SEXP, const matrix_info&);

// Functions to be called from R.

extern "C" {

SEXP colsum_subset(SEXP, SEXP);

SEXP colsum_exprs_subset(SEXP, SEXP, SEXP);

SEXP rowsum_exprs(SEXP, SEXP);

SEXP calc_top_features(SEXP, SEXP, SEXP);

SEXP negative_counts(SEXP);

SEXP missing_exprs(SEXP);

SEXP calc_exprs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif

