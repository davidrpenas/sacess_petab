#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<sunindextype, 9> dwdx_colptrs_Elowitz_Nature2000_ = {
    0, 3, 5, 7, 8, 10, 12, 14, 16
};

void dwdx_colptrs_Elowitz_Nature2000(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Elowitz_Nature2000_));
}
} // namespace model_Elowitz_Nature2000
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<sunindextype, 16> dwdx_rowvals_Elowitz_Nature2000_ = {
    8, 13, 15, 9, 14, 10, 12, 11, 0, 4, 1, 5, 2, 6, 3, 7
};

void dwdx_rowvals_Elowitz_Nature2000(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Elowitz_Nature2000_));
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
#include "dwdx.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void dwdx_Elowitz_Nature2000(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v9_v_8_dX_protein = 1.0*M_LN2/tau_prot;  // dwdx[0]
        dflux_v10_v_9_dY_protein = 1.0*M_LN2/tau_prot;  // dwdx[3]
        dflux_v11_v_10_dZ_protein = 1.0*M_LN2/tau_prot;  // dwdx[5]
        dflux_v12_v_11_dGFP = 1.0*M_LN2/tau_prot_GFP;  // dwdx[7]
        dflux_v1_v_0_dX_mRNA = 1.0*M_LN2/tau_mRNA;  // dwdx[8]
        dflux_v5_v_4_dX_mRNA = 1.0*eff/tau_mRNA;  // dwdx[9]
        dflux_v2_v_1_dY_mRNA = 1.0*M_LN2/tau_mRNA;  // dwdx[10]
        dflux_v6_v_5_dY_mRNA = 1.0*eff/tau_mRNA;  // dwdx[11]
        dflux_v3_v_2_dZ_mRNA = 1.0*M_LN2/tau_mRNA;  // dwdx[12]
        dflux_v7_v_6_dZ_mRNA = 1.0*eff/tau_mRNA;  // dwdx[13]
        dflux_v4_v_3_dGFP_mRNA = 1.0*M_LN2/tau_mRNA_GFP;  // dwdx[14]
        dflux_v8_v_7_dGFP_mRNA = 1.0*eff_GFP/tau_mRNA_GFP;  // dwdx[15]
    }

    // dynamic expressions
    dflux_v14_v_13_dX_protein = -std::pow(KM, n_Hill)*std::pow(X_protein, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill), 2);  // dwdx[1]
    dflux_v16_v_15_dX_protein = -std::pow(KM, n_Hill)*std::pow(X_protein, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(X_protein, n_Hill), 2);  // dwdx[2]
    dflux_v15_v_14_dY_protein = -std::pow(KM, n_Hill)*std::pow(Y_protein, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(Y_protein, n_Hill), 2);  // dwdx[4]
    dflux_v13_v_12_dZ_protein = -std::pow(KM, n_Hill)*std::pow(Z_protein, n_Hill - 1)*n_Hill*(60.0*tps_active - 60.0*tps_repr)/std::pow(std::pow(KM, n_Hill) + std::pow(Z_protein, n_Hill), 2);  // dwdx[6]
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
