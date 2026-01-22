#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void x_rdata_Armistead_CellDeathDis2024(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Cer;
    x_rdata[1] = Sphingo;
    x_rdata[2] = S1P;
    x_rdata[3] = Sphinga;
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
