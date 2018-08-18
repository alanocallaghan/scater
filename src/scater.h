#ifndef SCATER_H
#define SCATER_H

#include "Rcpp.h"
#include "beachmat/integer_matrix.h"
#include "beachmat/numeric_matrix.h"

#include <algorithm>
#include <deque>
#include <stdexcept>
#include <functional>
#include <cmath>

// Functions to be called from R.

extern "C" {

SEXP norm_exprs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP ave_exprs(SEXP, SEXP, SEXP, SEXP);


SEXP combined_qc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP top_cumprop(SEXP, SEXP);


SEXP sum_counts(SEXP, SEXP);


SEXP row_vars(SEXP, SEXP, SEXP);

SEXP col_vars(SEXP, SEXP, SEXP);

SEXP row_sums(SEXP, SEXP, SEXP);

SEXP col_sums(SEXP, SEXP, SEXP);

SEXP row_above(SEXP, SEXP, SEXP, SEXP);

SEXP col_above(SEXP, SEXP, SEXP, SEXP);

}

#include "utils.h"

#endif

