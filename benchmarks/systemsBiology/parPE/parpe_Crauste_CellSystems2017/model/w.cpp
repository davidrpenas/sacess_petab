#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void w_Crauste_CellSystems2017(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    flux_v1_v_0 = Naive*mu_N;  // w[0]
    flux_v2_v_1 = Naive*Pathogen*delta_NE;  // w[1]
    flux_v3_v_2 = EarlyEffector*Pathogen*rho_E;  // w[2]
    flux_v4_v_3 = std::pow(EarlyEffector, 2)*mu_EE;  // w[3]
    flux_v5_v_4 = EarlyEffector*delta_EL;  // w[4]
    flux_v6_v_5 = std::pow(LateEffector, 2)*mu_LL;  // w[5]
    flux_v7_v_6 = EarlyEffector*LateEffector*mu_LE;  // w[6]
    flux_v8_v_7 = LateEffector*delta_LM;  // w[7]
    flux_v9_v_8 = std::pow(Pathogen, 2)*rho_P;  // w[8]
    flux_v10_v_9 = EarlyEffector*Pathogen*mu_PE;  // w[9]
    flux_v11_v_10 = LateEffector*Pathogen*mu_PL;  // w[10]
    flux_v12_v_11 = Pathogen*mu_P;  // w[11]
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
