#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 11> dwdx_colptrs_Giordano_Nature2020_ = {
    0, 1, 6, 11, 16, 21, 25, 26, 27, 27, 27
};

void dwdx_colptrs_Giordano_Nature2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 27> dwdx_rowvals_Giordano_Nature2020_ = {
    17, 1, 17, 18, 19, 20, 0, 1, 17, 21, 22, 1, 17, 23, 24, 25, 0, 1, 17, 26, 27, 0, 1, 28, 29, 16, 16
};

void dwdx_rowvals_Giordano_Nature2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Giordano_Nature2020_));
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
#include "dwdx.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void dwdx_Giordano_Nature2020(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dCurrentTotalInfected_dInfected = 1;  // dwdx[1]
        dCurrentDiagnosedInfected_dDiagnosed = 1;  // dwdx[6]
        dCurrentTotalInfected_dDiagnosed = 1;  // dwdx[7]
        dCurrentTotalInfected_dAiling = 1;  // dwdx[11]
        dflux_Ailing_to_Recognised_dAiling = 1.0*theta;  // dwdx[13]
        dCurrentDiagnosedInfected_dRecognized = 1;  // dwdx[16]
        dCurrentTotalInfected_dRecognized = 1;  // dwdx[17]
        dCurrentDiagnosedInfected_dThreatened = 1;  // dwdx[21]
        dCurrentTotalInfected_dThreatened = 1;  // dwdx[22]
        dflux_Threatened_to_Extinct_dThreatened = 1.0*tau;  // dwdx[23]
        dCumulativeInfected_dHealed = 1;  // dwdx[25]
        dCumulativeInfected_dExtinct = 1;  // dwdx[26]
    }

    // dynamic expressions
    dflux_Susceptible_to_Infected_dSusceptible = 1.0*Ailing*gamma + 1.0*Diagnosed*beta + 1.0*Infected*alpha + 1.0*Recognized*delta;  // dwdx[0]
    dflux_Susceptible_to_Infected_dInfected = 1.0*Susceptible*alpha;  // dwdx[2]
    dflux_Infected_to_Diagnosed_dInfected = 1.0*epsilon;  // dwdx[3]
    dflux_Infected_to_Ailing_dInfected = 1.0*zeta;  // dwdx[4]
    dflux_Infected_to_Healed_dInfected = 1.0*lam;  // dwdx[5]
    dflux_Susceptible_to_Infected_dDiagnosed = 1.0*Susceptible*beta;  // dwdx[8]
    dflux_Diagnosed_to_Recognized_dDiagnosed = 1.0*eta;  // dwdx[9]
    dflux_Diagnosed_to_Healed_dDiagnosed = 1.0*rho;  // dwdx[10]
    dflux_Susceptible_to_Infected_dAiling = 1.0*Susceptible*gamma;  // dwdx[12]
    dflux_Ailing_to_Healed_dAiling = 1.0*kappa;  // dwdx[14]
    dflux_Ailing_to_Threatened_dAiling = 1.0*mu;  // dwdx[15]
    dflux_Susceptible_to_Infected_dRecognized = 1.0*Susceptible*delta;  // dwdx[18]
    dflux_Recognised_to_Threatened_dRecognized = 1.0*nu;  // dwdx[19]
    dflux_Recognised_to_Healed_dRecognized = 1.0*xi;  // dwdx[20]
    dflux_Threatened_to_Healed_dThreatened = 1.0*sigma;  // dwdx[24]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
