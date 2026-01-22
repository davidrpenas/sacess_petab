#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<sunindextype, 22> dwdp_colptrs_Elowitz_Nature2000_ = {
    0, 4, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 12, 18, 20, 23, 24, 28, 32, 32, 32, 32
};

void dwdp_colptrs_Elowitz_Nature2000(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Elowitz_Nature2000_));
}
} // namespace model_Elowitz_Nature2000
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<sunindextype, 32> dwdp_rowvals_Elowitz_Nature2000_ = {
    12, 13, 14, 15, 4, 5, 6, 7, 12, 13, 14, 15, 0, 1, 2, 4, 5, 6, 3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 12, 13, 14, 15
};

void dwdp_rowvals_Elowitz_Nature2000(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Elowitz_Nature2000_));
}
} // namespace model_Elowitz_Nature2000
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "dwdp.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void dwdp_Elowitz_Nature2000(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v13_v_12_dKM = std::pow(KM, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill)) - std::pow(KM, 2*n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill), 2);  // dwdp[0]
    dflux_v14_v_13_dKM = std::pow(KM, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill)) - std::pow(KM, 2*n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill), 2);  // dwdp[1]
    dflux_v15_v_14_dKM = std::pow(KM, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill)) - std::pow(KM, 2*n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill), 2);  // dwdp[2]
    dflux_v16_v_15_dKM = std::pow(KM, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill)) - std::pow(KM, 2*n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill), 2);  // dwdp[3]
    dflux_v5_v_4_deff = 1.0*X_mRNA/tau_mRNA;  // dwdp[4]
    dflux_v6_v_5_deff = 1.0*Y_mRNA/tau_mRNA;  // dwdp[5]
    dflux_v7_v_6_deff = 1.0*Z_mRNA/tau_mRNA;  // dwdp[6]
    dflux_v8_v_7_deff_GFP = 1.0*GFP_mRNA/tau_mRNA_GFP;  // dwdp[7]
    dflux_v13_v_12_dn_Hill = -std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)*(((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    )) + ((Z_protein == 0) ? (
       Z_protein
    )
    : (
       std::pow(Z_protein, n_Hill)*std::log(Z_protein)
    )))/std::pow(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill), 2) + (60.0*tps_active - 60.0*tps_repr)*((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    ))/(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill));  // dwdp[8]
    dflux_v14_v_13_dn_Hill = -std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)*(((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    )) + ((X_protein == 0) ? (
       X_protein
    )
    : (
       std::pow(X_protein, n_Hill)*std::log(X_protein)
    )))/std::pow(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill), 2) + (60.0*tps_active - 60.0*tps_repr)*((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    ))/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill));  // dwdp[9]
    dflux_v15_v_14_dn_Hill = -std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)*(((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    )) + ((Y_protein == 0) ? (
       Y_protein
    )
    : (
       std::pow(Y_protein, n_Hill)*std::log(Y_protein)
    )))/std::pow(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill), 2) + (60.0*tps_active - 60.0*tps_repr)*((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    ))/(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill));  // dwdp[10]
    dflux_v16_v_15_dn_Hill = -std::pow(KM, n_Hill)*(60.0*tps_active - 60.0*tps_repr)*(((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    )) + ((X_protein == 0) ? (
       X_protein
    )
    : (
       std::pow(X_protein, n_Hill)*std::log(X_protein)
    )))/std::pow(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill), 2) + (60.0*tps_active - 60.0*tps_repr)*((KM == 0) ? (
       KM
    )
    : (
       std::pow(KM, n_Hill)*std::log(KM)
    ))/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill));  // dwdp[11]
    dflux_v1_v_0_dtau_mRNA = -1.0*X_mRNA*M_LN2/std::pow(tau_mRNA, 2);  // dwdp[12]
    dflux_v2_v_1_dtau_mRNA = -1.0*Y_mRNA*M_LN2/std::pow(tau_mRNA, 2);  // dwdp[13]
    dflux_v3_v_2_dtau_mRNA = -1.0*Z_mRNA*M_LN2/std::pow(tau_mRNA, 2);  // dwdp[14]
    dflux_v5_v_4_dtau_mRNA = -1.0*X_mRNA*eff/std::pow(tau_mRNA, 2);  // dwdp[15]
    dflux_v6_v_5_dtau_mRNA = -1.0*Y_mRNA*eff/std::pow(tau_mRNA, 2);  // dwdp[16]
    dflux_v7_v_6_dtau_mRNA = -1.0*Z_mRNA*eff/std::pow(tau_mRNA, 2);  // dwdp[17]
    dflux_v4_v_3_dtau_mRNA_GFP = -1.0*GFP_mRNA*M_LN2/std::pow(tau_mRNA_GFP, 2);  // dwdp[18]
    dflux_v8_v_7_dtau_mRNA_GFP = -1.0*GFP_mRNA*eff_GFP/std::pow(tau_mRNA_GFP, 2);  // dwdp[19]
    dflux_v9_v_8_dtau_prot = -1.0*X_protein*M_LN2/std::pow(tau_prot, 2);  // dwdp[20]
    dflux_v10_v_9_dtau_prot = -1.0*Y_protein*M_LN2/std::pow(tau_prot, 2);  // dwdp[21]
    dflux_v11_v_10_dtau_prot = -1.0*Z_protein*M_LN2/std::pow(tau_prot, 2);  // dwdp[22]
    dflux_v12_v_11_dtau_prot_GFP = -1.0*GFP*M_LN2/std::pow(tau_prot_GFP, 2);  // dwdp[23]
    dflux_v13_v_12_dtps_active = 60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill));  // dwdp[24]
    dflux_v14_v_13_dtps_active = 60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill));  // dwdp[25]
    dflux_v15_v_14_dtps_active = 60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill));  // dwdp[26]
    dflux_v16_v_15_dtps_active = 60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill));  // dwdp[27]
    dflux_v13_v_12_dtps_repr = -60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill)) + 60.0;  // dwdp[28]
    dflux_v14_v_13_dtps_repr = -60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill)) + 60.0;  // dwdp[29]
    dflux_v15_v_14_dtps_repr = -60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill)) + 60.0;  // dwdp[30]
    dflux_v16_v_15_dtps_repr = -60.0*std::pow(KM, n_Hill)/(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill)) + 60.0;  // dwdp[31]
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
