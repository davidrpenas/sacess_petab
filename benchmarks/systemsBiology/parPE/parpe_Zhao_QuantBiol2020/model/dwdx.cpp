#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 5> dwdx_colptrs_Zhao_QuantBiol2020_ = {
    0, 1, 4, 6, 7
};

void dwdx_colptrs_Zhao_QuantBiol2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Zhao_QuantBiol2020_));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 7> dwdx_rowvals_Zhao_QuantBiol2020_ = {
    7, 0, 7, 8, 0, 9, 0
};

void dwdx_rowvals_Zhao_QuantBiol2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Zhao_QuantBiol2020_));
}
} // namespace model_Zhao_QuantBiol2020
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
namespace model_Zhao_QuantBiol2020 {

void dwdx_Zhao_QuantBiol2020(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dCumulative_Infected_dUnquarantined_Infected = 1;  // dwdx[1]
        dflux_Unquarantined_to_Quarantined_dUnquarantined_Infected = 1.0*gamma_1;  // dwdx[3]
        dCumulative_Infected_dQuarantined_Infected = 1;  // dwdx[4]
        dflux_Quarantined_to_Confirmed_dQuarantined_Infected = 1.0*beta;  // dwdx[5]
        dCumulative_Infected_dConfirmed_Infected = 1;  // dwdx[6]
    }

    // dynamic expressions
    dflux_Susceptible_to_Unquarantined_dSusceptible = 1.0*Unquarantined_Infected*alpha/Total_Pop;  // dwdx[0]
    dflux_Susceptible_to_Unquarantined_dUnquarantined_Infected = 1.0*Susceptible*alpha/Total_Pop;  // dwdx[2]
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
