#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void x_rdata_Borghans_BiophysChem1997(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Z_state;
    x_rdata[1] = Y_state;
    x_rdata[2] = A_state;
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
