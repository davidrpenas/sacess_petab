#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 11> dwdw_colptrs_Zhao_QuantBiol2020_ = {
    0, 0, 1, 2, 4, 5, 6, 7, 7, 7, 7
};

void dwdw_colptrs_Zhao_QuantBiol2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Zhao_QuantBiol2020_));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 7> dwdw_rowvals_Zhao_QuantBiol2020_ = {
    5, 7, 5, 8, 6, 7, 9
};

void dwdw_rowvals_Zhao_QuantBiol2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Zhao_QuantBiol2020_));
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
#include "dwdw.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void dwdw_Zhao_QuantBiol2020(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dalpha_dR = gamma_1;  // dwdw[0]
        dalpha_dgamma_1 = R;  // dwdw[2]
        dbeta_dgamma_2 = 1 - sigma;  // dwdw[4]
    }

    // dynamic expressions
    dflux_Susceptible_to_Unquarantined_dTotal_Pop = -1.0*Susceptible*Unquarantined_Infected*alpha/std::pow(Total_Pop, 2);  // dwdw[1]
    dflux_Unquarantined_to_Quarantined_dgamma_1 = 1.0*Unquarantined_Infected;  // dwdw[3]
    dflux_Susceptible_to_Unquarantined_dalpha = 1.0*Susceptible*Unquarantined_Infected/Total_Pop;  // dwdw[5]
    dflux_Quarantined_to_Confirmed_dbeta = 1.0*Quarantined_Infected;  // dwdw[6]
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
