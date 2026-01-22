#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void x_rdata_Zheng_PNAS2012(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = K27me0K36me0;
    x_rdata[1] = K27me0K36me1;
    x_rdata[2] = K27me1K36me0;
    x_rdata[3] = K27me0K36me2;
    x_rdata[4] = K27me1K36me1;
    x_rdata[5] = K27me2K36me0;
    x_rdata[6] = K27me0K36me3;
    x_rdata[7] = K27me1K36me2;
    x_rdata[8] = K27me2K36me1;
    x_rdata[9] = K27me3K36me0;
    x_rdata[10] = K27me1K36me3;
    x_rdata[11] = K27me2K36me2;
    x_rdata[12] = K27me3K36me1;
    x_rdata[13] = K27me2K36me3;
    x_rdata[14] = K27me3K36me2;
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
