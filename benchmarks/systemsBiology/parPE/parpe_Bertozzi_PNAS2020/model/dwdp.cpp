#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 6> dwdp_colptrs_Bertozzi_PNAS2020_ = {
    0, 2, 2, 3, 3, 3
};

void dwdp_colptrs_Bertozzi_PNAS2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Bertozzi_PNAS2020_));
}
} // namespace model_Bertozzi_PNAS2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 3> dwdp_rowvals_Bertozzi_PNAS2020_ = {
    0, 2, 0
};

void dwdp_rowvals_Bertozzi_PNAS2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Bertozzi_PNAS2020_));
}
} // namespace model_Bertozzi_PNAS2020
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
#include "dwdp.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void dwdp_Bertozzi_PNAS2020(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dbeta_N_dgamma_ = R0_/N_;  // dwdp[0]
        dbeta_N_dR0_ = gamma_/N_;  // dwdp[2]
    }

    // dynamic expressions
    dflux_recovery_dgamma_ = I_;  // dwdp[1]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
