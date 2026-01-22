#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void y_Zheng_PNAS2012(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = K27me0K36me0/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[1] = K27me0K36me1/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[2] = K27me0K36me2/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[3] = K27me0K36me3/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[4] = K27me1K36me0/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[5] = K27me1K36me1/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[6] = K27me1K36me2/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[7] = K27me1K36me3/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[8] = K27me2K36me0/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[9] = K27me2K36me1/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[10] = K27me2K36me2/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[11] = K27me2K36me3/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[12] = K27me3K36me0/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[13] = K27me3K36me1/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
    y[14] = K27me3K36me2/(K27me0K36me0 + K27me0K36me1 + K27me0K36me2 + K27me0K36me3 + K27me1K36me0 + K27me1K36me1 + K27me1K36me2 + K27me1K36me3 + K27me2K36me0 + K27me2K36me1 + K27me2K36me2 + K27me2K36me3 + K27me3K36me0 + K27me3K36me1 + K27me3K36me2);
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
