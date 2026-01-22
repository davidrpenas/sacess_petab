#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 10> dxdotdx_explicit_colptrs_Oliveira_NatCommun2021_ = {
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

void dxdotdx_explicit_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexptrs(gsl::make_span(dxdotdx_explicit_colptrs_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 1> dxdotdx_explicit_rowvals_Oliveira_NatCommun2021_ = {
    8
};

void dxdotdx_explicit_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &dxdotdx_explicit){
    dxdotdx_explicit.set_indexvals(gsl::make_span(dxdotdx_explicit_rowvals_Oliveira_NatCommun2021_));
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
#include "dxdotdx_explicit.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void dxdotdx_explicit_Oliveira_NatCommun2021(realtype *dxdotdx_explicit, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddCumulative_casesdt_dExposed = 1.0*exposed_to_symptomatic;  // dxdotdx_explicit[0]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
