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

SEXP calc_top_features(SEXP, SEXP, SEXP);

SEXP margin_summary(SEXP, SEXP, SEXP, SEXP);

SEXP calc_exprs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP calc_variance(SEXP, SEXP, SEXP);

SEXP downsample_matrix(SEXP, SEXP);

}

#include "template_methods.h"

#endif

