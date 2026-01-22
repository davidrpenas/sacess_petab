#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 13> dxdotdw_colptrs_Okuonghae_ChaosSolitonsFractals2020_ = {
    0, 0, 0, 2, 4, 6, 9, 12, 14, 16, 18, 20, 22
};

void dxdotdw_colptrs_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 22> dxdotdw_rowvals_Okuonghae_ChaosSolitonsFractals2020_ = {
    0, 1, 1, 3, 1, 2, 3, 4, 7, 2, 4, 7, 3, 5, 2, 5, 4, 5, 3, 6, 4, 6
};

void dxdotdw_rowvals_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Okuonghae_ChaosSolitonsFractals2020_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void dxdotdw_Okuonghae_ChaosSolitonsFractals2020(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddsusceptibledt_dflux_susceptible_to_exposed = -1;  // dxdotdw[0]
    ddexposeddt_dflux_susceptible_to_exposed = 1;  // dxdotdw[1]
    ddexposeddt_dflux_exposed_to_symptomatic = -1;  // dxdotdw[2]
    ddsymptomaticdt_dflux_exposed_to_symptomatic = 1;  // dxdotdw[3]
    ddexposeddt_dflux_exposed_to_asymptopmatic = -1;  // dxdotdw[4]
    ddasymptomaticdt_dflux_exposed_to_asymptopmatic = 1;  // dxdotdw[5]
    ddsymptomaticdt_dflux_symptomatic_to_detected = -1;  // dxdotdw[6]
    dddetecteddt_dflux_symptomatic_to_detected = 1;  // dxdotdw[7]
    dddetected_cumulativedt_dflux_symptomatic_to_detected = 1;  // dxdotdw[8]
    ddasymptomaticdt_dflux_asymptomatic_to_detected = -1;  // dxdotdw[9]
    dddetecteddt_dflux_asymptomatic_to_detected = 1;  // dxdotdw[10]
    dddetected_cumulativedt_dflux_asymptomatic_to_detected = 1;  // dxdotdw[11]
    ddsymptomaticdt_dflux_symptomatic_to_recovered = -1;  // dxdotdw[12]
    ddrecovereddt_dflux_symptomatic_to_recovered = 1;  // dxdotdw[13]
    ddasymptomaticdt_dflux_asymptomatic_to_recovered = -1;  // dxdotdw[14]
    ddrecovereddt_dflux_asymptomatic_to_recovered = 1;  // dxdotdw[15]
    dddetecteddt_dflux_detected_to_recovered = -1;  // dxdotdw[16]
    ddrecovereddt_dflux_detected_to_recovered = 1;  // dxdotdw[17]
    ddsymptomaticdt_dflux_symptomatic_to_deceased = -1;  // dxdotdw[18]
    dddeceaseddt_dflux_symptomatic_to_deceased = 1;  // dxdotdw[19]
    dddetecteddt_dflux_detected_to_deceased = -1;  // dxdotdw[20]
    dddeceaseddt_dflux_detected_to_deceased = 1;  // dxdotdw[21]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
