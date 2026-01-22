#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 23> dxdotdw_colptrs_Oliveira_NatCommun2021_ = {
    0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25
};

void dxdotdw_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 25> dxdotdw_rowvals_Oliveira_NatCommun2021_ = {
    8, 0, 1, 0, 2, 0, 3, 2, 4, 3, 4, 3, 5, 4, 5, 3, 6, 5, 6, 5, 6, 5, 7, 6, 7
};

void dxdotdw_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
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
#include "dxdotdw.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void dxdotdw_Oliveira_NatCommun2021(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddCumulative_casesdt_dexposed_to_symptomatic = 1.0*Exposed;  // dxdotdw[0]
    ddExposeddt_dflux_Susceptible_to_Exposed = 1;  // dxdotdw[1]
    ddSusceptibledt_dflux_Susceptible_to_Exposed = -1;  // dxdotdw[2]
    ddExposeddt_dflux_Exposed_to_Asymptomatic = -1;  // dxdotdw[3]
    ddAsymptomaticdt_dflux_Exposed_to_Asymptomatic = 1;  // dxdotdw[4]
    ddExposeddt_dflux_Exposed_to_Symptomatic = -1;  // dxdotdw[5]
    ddSymptomaticdt_dflux_Exposed_to_Symptomatic = 1;  // dxdotdw[6]
    ddAsymptomaticdt_dflux_Asymptomatic_to_Recovered = -1;  // dxdotdw[7]
    ddRecovereddt_dflux_Asymptomatic_to_Recovered = 1;  // dxdotdw[8]
    ddSymptomaticdt_dflux_Symptomatic_to_Recovered = -1;  // dxdotdw[9]
    ddRecovereddt_dflux_Symptomatic_to_Recovered = 1;  // dxdotdw[10]
    ddSymptomaticdt_dflux_Symptomatic_to_Hospitalized = -1;  // dxdotdw[11]
    ddHospitaldt_dflux_Symptomatic_to_Hospitalized = 1;  // dxdotdw[12]
    ddRecovereddt_dflux_Hospitalized_to_Recovered = 1;  // dxdotdw[13]
    ddHospitaldt_dflux_Hospitalized_to_Recovered = -1;  // dxdotdw[14]
    ddSymptomaticdt_dflux_Symptomatic_to_ICU = -1;  // dxdotdw[15]
    ddICUdt_dflux_Symptomatic_to_ICU = 1;  // dxdotdw[16]
    ddHospitaldt_dflux_Hospitalized_to_ICU = -1;  // dxdotdw[17]
    ddICUdt_dflux_Hospitalized_to_ICU = 1;  // dxdotdw[18]
    ddHospitaldt_dflux_ICU_to_Hospitalized = 1;  // dxdotdw[19]
    ddICUdt_dflux_ICU_to_Hospitalized = -1;  // dxdotdw[20]
    ddHospitaldt_dflux_Hospitalized_to_Deaths = -1;  // dxdotdw[21]
    ddDeathsdt_dflux_Hospitalized_to_Deaths = 1;  // dxdotdw[22]
    ddICUdt_dflux_ICU_to_Deaths = -1;  // dxdotdw[23]
    ddDeathsdt_dflux_ICU_to_Deaths = 1;  // dxdotdw[24]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
