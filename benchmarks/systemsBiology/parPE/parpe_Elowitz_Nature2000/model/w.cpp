#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void w_Elowitz_Nature2000(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    flux_v1_v_0 = 1.0*X_mRNA*M_LN2/tau_mRNA;  // w[0]
    flux_v2_v_1 = 1.0*Y_mRNA*M_LN2/tau_mRNA;  // w[1]
    flux_v3_v_2 = 1.0*Z_mRNA*M_LN2/tau_mRNA;  // w[2]
    flux_v4_v_3 = 1.0*GFP_mRNA*M_LN2/tau_mRNA_GFP;  // w[3]
    flux_v5_v_4 = 1.0*X_mRNA*eff/tau_mRNA;  // w[4]
    flux_v6_v_5 = 1.0*Y_mRNA*eff/tau_mRNA;  // w[5]
    flux_v7_v_6 = 1.0*Z_mRNA*eff/tau_mRNA;  // w[6]
    flux_v8_v_7 = 1.0*GFP_mRNA*eff_GFP/tau_mRNA_GFP;  // w[7]
    flux_v9_v_8 = 1.0*X_protein*M_LN2/tau_prot;  // w[8]
    flux_v10_v_9 = 1.0*Y_protein*M_LN2/tau_prot;  // w[9]
    flux_v11_v_10 = 1.0*Z_protein*M_LN2/tau_prot;  // w[10]
    flux_v12_v_11 = 1.0*GFP*M_LN2/tau_prot_GFP;  // w[11]
    flux_v13_v_12 = std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill)) + 60.0*tps_repr;  // w[12]
    flux_v14_v_13 = std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill)) + 60.0*tps_repr;  // w[13]
    flux_v15_v_14 = std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill)) + 60.0*tps_repr;  // w[14]
    flux_v16_v_15 = std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill)) + 60.0*tps_repr;  // w[15]
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
