#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 23> dwdw_colptrs_Oliveira_NatCommun2021_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
};

void dwdw_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 10> dwdw_rowvals_Oliveira_NatCommun2021_ = {
    10, 11, 12, 20, 18, 16, 21, 15, 17, 14
};

void dwdw_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Oliveira_NatCommun2021_));
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
#include "dwdw.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void dwdw_Oliveira_NatCommun2021(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_Susceptible_to_Exposed_dbeta = 1.0*Susceptible*(Asymptomatic*delta_ + Symptomatic)/population;  // dwdw[0]
    dflux_Exposed_to_Asymptomatic_dexposed_to_asymptomatic = 1.0*Exposed;  // dwdw[1]
    dflux_Exposed_to_Symptomatic_dexposed_to_symptomatic = 1.0*Exposed;  // dwdw[2]
    dflux_Hospitalized_to_Deaths_dhospitalized_to_deaths = 1.0*Hospital;  // dwdw[3]
    dflux_Hospitalized_to_ICU_dhospitalized_to_icu = 1.0*Hospital;  // dwdw[4]
    dflux_Hospitalized_to_Recovered_dhospitalized_to_recovered = 1.0*Hospital;  // dwdw[5]
    dflux_ICU_to_Deaths_dicu_to_deaths = 1.0*ICU;  // dwdw[6]
    dflux_Symptomatic_to_Hospitalized_dsymptomatic_to_hospitalized = 1.0*Symptomatic;  // dwdw[7]
    dflux_Symptomatic_to_ICU_dsymptomatic_to_icu = 1.0*Symptomatic;  // dwdw[8]
    dflux_Symptomatic_to_Recovered_dsymptomatic_to_recovered = 1.0*Symptomatic;  // dwdw[9]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
