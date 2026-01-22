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
namespace model_Laske_PLOSComputBiol2019 {

void y_Laske_PLOSComputBiol2019(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = R_M5;
    y[1] = R_V_seg_tot + observableParameter1_RVSegTot;
    y[2] = R_C_seg_tot;
    y[3] = V_rel;
    y[4] = observableParameter1_IntNucOffset;
    y[5] = F_rnp_nuc + observableParameter1_FracNucInt_1;
    y[6] = F_rnp_nuc + observableParameter1_FracNucInt_2;
    y[7] = F_rnp_nuc + observableParameter1_FracNucInt_3;
    y[8] = F_rnp_nuc + observableParameter1_FracNucInt_4;
    y[9] = F_rnp_nuc + observableParameter1_FracNucInt_5;
    y[10] = F_rnp_nuc + observableParameter1_FracNucInt_6;
    y[11] = F_rnp_nuc + observableParameter1_FracNucInt_7;
    y[12] = F_rnp_nuc + observableParameter1_FracNucInt_8;
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
