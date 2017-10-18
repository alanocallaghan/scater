#include "scater.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define REGISTER(x, i) {#x, (DL_FUNC) &x, i}

extern "C" {

static const R_CallMethodDef all_call_entries[] = {
    REGISTER(calc_exprs, 7),
    REGISTER(margin_summary, 4),
    REGISTER(calc_top_features, 3),
    REGISTER(calc_variance, 3),
    {NULL, NULL, 0}
};

void attribute_visible R_init_scater(DllInfo *dll) {
    R_registerRoutines(dll, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

}
