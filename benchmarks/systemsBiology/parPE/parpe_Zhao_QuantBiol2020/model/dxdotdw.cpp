#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 11> dxdotdw_colptrs_Zhao_QuantBiol2020_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 6
};

void dxdotdw_colptrs_Zhao_QuantBiol2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Zhao_QuantBiol2020_));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 6> dxdotdw_rowvals_Zhao_QuantBiol2020_ = {
    0, 1, 1, 2, 2, 3
};

void dxdotdw_rowvals_Zhao_QuantBiol2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Zhao_QuantBiol2020_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void dxdotdw_Zhao_QuantBiol2020(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSusceptibledt_dflux_Susceptible_to_Unquarantined = -1;  // dxdotdw[0]
    ddUnquarantined_Infecteddt_dflux_Susceptible_to_Unquarantined = 1;  // dxdotdw[1]
    ddUnquarantined_Infecteddt_dflux_Unquarantined_to_Quarantined = -1;  // dxdotdw[2]
    ddQuarantined_Infecteddt_dflux_Unquarantined_to_Quarantined = 1;  // dxdotdw[3]
    ddQuarantined_Infecteddt_dflux_Quarantined_to_Confirmed = -1;  // dxdotdw[4]
    ddConfirmed_Infecteddt_dflux_Quarantined_to_Confirmed = 1;  // dxdotdw[5]
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
