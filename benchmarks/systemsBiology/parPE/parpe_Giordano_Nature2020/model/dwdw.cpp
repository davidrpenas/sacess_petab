#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 31> dwdw_colptrs_Giordano_Nature2020_ = {
    0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
};

void dwdw_colptrs_Giordano_Nature2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 15> dwdw_rowvals_Giordano_Nature2020_ = {
    16, 17, 17, 17, 18, 21, 17, 24, 20, 25, 26, 22, 29, 27, 19
};

void dwdw_rowvals_Giordano_Nature2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dwdw.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void dwdw_Giordano_Nature2020(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dCumulativeInfected_dCurrentTotalInfected = 1;  // dwdw[0]
    }

    // dynamic expressions
    dflux_Susceptible_to_Infected_dalpha = 1.0*Infected*Susceptible;  // dwdw[1]
    dflux_Susceptible_to_Infected_dbeta = 1.0*Diagnosed*Susceptible;  // dwdw[2]
    dflux_Susceptible_to_Infected_ddelta = 1.0*Recognized*Susceptible;  // dwdw[3]
    dflux_Infected_to_Diagnosed_depsilon = 1.0*Infected;  // dwdw[4]
    dflux_Diagnosed_to_Recognized_deta = 1.0*Diagnosed;  // dwdw[5]
    dflux_Susceptible_to_Infected_dgamma = 1.0*Ailing*Susceptible;  // dwdw[6]
    dflux_Ailing_to_Healed_dkappa = 1.0*Ailing;  // dwdw[7]
    dflux_Infected_to_Healed_dlam = 1.0*Infected;  // dwdw[8]
    dflux_Ailing_to_Threatened_dmu = 1.0*Ailing;  // dwdw[9]
    dflux_Recognised_to_Threatened_dnu = 1.0*Recognized;  // dwdw[10]
    dflux_Diagnosed_to_Healed_drho = 1.0*Diagnosed;  // dwdw[11]
    dflux_Threatened_to_Healed_dsigma = 1.0*Threatened;  // dwdw[12]
    dflux_Recognised_to_Healed_dxi = 1.0*Recognized;  // dwdw[13]
    dflux_Infected_to_Ailing_dzeta = 1.0*Infected;  // dwdw[14]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
