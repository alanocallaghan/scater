#include "scater.h"
#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(calc_exprs, 6),
    REGISTER(colsum_subset, 2),
    REGISTER(colsum_exprs_subset, 3),
    REGISTER(rowsum_exprs, 2),
    REGISTER(negative_counts, 1),
    REGISTER(missing_exprs, 1),
    REGISTER(calc_top_features, 3),
    {NULL, NULL, 0}
};

void attribute_visible R_init_scater(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
