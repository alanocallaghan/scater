#ifndef SCATER_H
#define SCATER_H

#include "Rcpp.h"

// Functions to be called from R.

extern "C" {

SEXP norm_exprs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP ave_exprs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


SEXP combined_qc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP top_cumprop(SEXP, SEXP);


SEXP sum_counts(SEXP, SEXP);


SEXP row_above(SEXP, SEXP, SEXP, SEXP);

SEXP col_above(SEXP, SEXP, SEXP, SEXP);

}

#endif

