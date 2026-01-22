#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 4> dwdw_colptrs_Bertozzi_PNAS2020_ = {
    0, 1, 1, 1
};

void dwdw_colptrs_Bertozzi_PNAS2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Bertozzi_PNAS2020_));
}
} // namespace model_Bertozzi_PNAS2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<sunindextype, 1> dwdw_rowvals_Bertozzi_PNAS2020_ = {
    1
};

void dwdw_rowvals_Bertozzi_PNAS2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Bertozzi_PNAS2020_));
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
#include "dwdw.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void dwdw_Bertozzi_PNAS2020(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_infection_dbeta_N = I_*S_;  // dwdw[0]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
