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
namespace model_Armistead_CellDeathDis2024 {

void w_Armistead_CellDeathDis2024(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        flux_R1 = k00*(S_on*alpha_cer + 1);  // w[0]
    }

    // dynamic expressions
    flux_R2 = Sphinga*k0;  // w[1]
    flux_R3 = Cer*k1;  // w[2]
    flux_R4 = Sphingo*k2;  // w[3]
    flux_R5 = Cer*k_d;  // w[4]
    flux_R6 = Sphingo*k3*(-S_on*alpha_hai1a + 1);  // w[5]
    flux_R7 = S1P*k4;  // w[6]
    flux_R8 = S1P*k5;  // w[7]
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
