#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void x_rdata_Crauste_CellSystems2017(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Naive;
    x_rdata[1] = EarlyEffector;
    x_rdata[2] = LateEffector;
    x_rdata[3] = Memory;
    x_rdata[4] = Pathogen;
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
