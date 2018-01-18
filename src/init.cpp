#include "scater.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(calc_exprs, 7),
    REGISTER(calc_top_features, 3),
    REGISTER(row_vars, 3),
    REGISTER(col_vars, 3),
    REGISTER(row_sums, 3),
    REGISTER(col_sums, 3),
    REGISTER(row_above, 4),
    REGISTER(col_above, 4),
    {NULL, NULL, 0}
};

void attribute_visible R_init_scater(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
