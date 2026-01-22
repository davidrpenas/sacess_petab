#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 4> dxdotdw_colptrs_Bertozzi_PNAS2020_ = {
    0, 0, 2, 4
};

void dxdotdw_colptrs_Bertozzi_PNAS2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Bertozzi_PNAS2020_));
}
} // namespace model_Bertozzi_PNAS2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 4> dxdotdw_rowvals_Bertozzi_PNAS2020_ = {
    0, 2, 0, 1
};

void dxdotdw_rowvals_Bertozzi_PNAS2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Bertozzi_PNAS2020_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void dxdotdw_Bertozzi_PNAS2020(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddI_dt_dflux_infection = 1;  // dxdotdw[0]
    ddS_dt_dflux_infection = -1;  // dxdotdw[1]
    ddI_dt_dflux_recovery = -1;  // dxdotdw[2]
    ddR_dt_dflux_recovery = 1;  // dxdotdw[3]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
