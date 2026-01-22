#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void x_solver_Zheng_PNAS2012(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = K27me0K36me0;
    x_solver[1] = K27me0K36me1;
    x_solver[2] = K27me1K36me0;
    x_solver[3] = K27me0K36me2;
    x_solver[4] = K27me1K36me1;
    x_solver[5] = K27me2K36me0;
    x_solver[6] = K27me0K36me3;
    x_solver[7] = K27me1K36me2;
    x_solver[8] = K27me2K36me1;
    x_solver[9] = K27me3K36me0;
    x_solver[10] = K27me1K36me3;
    x_solver[11] = K27me2K36me2;
    x_solver[12] = K27me3K36me1;
    x_solver[13] = K27me2K36me3;
    x_solver[14] = K27me3K36me2;
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
