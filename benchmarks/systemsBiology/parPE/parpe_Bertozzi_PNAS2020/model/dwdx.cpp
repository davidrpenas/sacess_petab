#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 4> dwdx_colptrs_Bertozzi_PNAS2020_ = {
    0, 2, 2, 3
};

void dwdx_colptrs_Bertozzi_PNAS2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Bertozzi_PNAS2020_));
}
} // namespace model_Bertozzi_PNAS2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 3> dwdx_rowvals_Bertozzi_PNAS2020_ = {
    1, 2, 1
};

void dwdx_rowvals_Bertozzi_PNAS2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Bertozzi_PNAS2020_));
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
#include "dwdx.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void dwdx_Bertozzi_PNAS2020(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_recovery_dI_ = gamma_;  // dwdx[1]
    }

    // dynamic expressions
    dflux_infection_dI_ = S_*beta_N;  // dwdx[0]
    dflux_infection_dS_ = I_*beta_N;  // dwdx[2]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
