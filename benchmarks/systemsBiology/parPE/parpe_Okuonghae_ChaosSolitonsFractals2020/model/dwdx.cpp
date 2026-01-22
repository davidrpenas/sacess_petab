#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 9> dwdx_colptrs_Okuonghae_ChaosSolitonsFractals2020_ = {
    0, 2, 5, 9, 14, 18, 19, 19, 19
};

void dwdx_colptrs_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 19> dwdx_rowvals_Okuonghae_ChaosSolitonsFractals2020_ = {
    0, 2, 0, 3, 4, 0, 1, 6, 8, 0, 1, 5, 7, 10, 0, 1, 9, 11, 0
};

void dwdx_rowvals_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void dwdx_Okuonghae_ChaosSolitonsFractals2020(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dtotal_pop_dsusceptible = 1;  // dwdx[0]
        dtotal_pop_dexposed = 1;  // dwdx[2]
        dflux_exposed_to_symptomatic_dexposed = 1.0*sigma*(1 - nu);  // dwdx[3]
        dflux_exposed_to_asymptopmatic_dexposed = 1.0*nu*sigma;  // dwdx[4]
        dtotal_pop_dasymptomatic = 1;  // dwdx[5]
        dflux_asymptomatic_to_detected_dasymptomatic = 1.0*theta;  // dwdx[7]
        dflux_asymptomatic_to_recovered_dasymptomatic = 1.0*gamma_a;  // dwdx[8]
        dtotal_pop_dsymptomatic = 1;  // dwdx[9]
        dflux_symptomatic_to_detected_dsymptomatic = 1.0*psi;  // dwdx[11]
        dflux_symptomatic_to_recovered_dsymptomatic = 1.0*gamma_0;  // dwdx[12]
        dflux_symptomatic_to_deceased_dsymptomatic = 1.0*d_0;  // dwdx[13]
        dtotal_pop_ddetected = 1;  // dwdx[14]
        dflux_detected_to_recovered_ddetected = 1.0*gamma_i;  // dwdx[16]
        dflux_detected_to_deceased_ddetected = 1.0*d_D;  // dwdx[17]
        dtotal_pop_drecovered = 1;  // dwdx[18]
    }

    // dynamic expressions
    dflux_susceptible_to_exposed_dsusceptible = 1.0*lam;  // dwdx[1]
    dlam_dasymptomatic = alpha*transmission_rate_effective*(1 - delta)*(1 - eps)/(-detected + total_pop);  // dwdx[6]
    dlam_dsymptomatic = transmission_rate_effective*(1 - delta)*(1 - eps)/(-detected + total_pop);  // dwdx[10]
    dlam_ddetected = transmission_rate_effective*(1 - delta)*(1 - eps)*(alpha*asymptomatic + symptomatic)/std::pow(-detected + total_pop, 2);  // dwdx[15]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
