#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 13> dwdw_colptrs_Okuonghae_ChaosSolitonsFractals2020_ = {
    0, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2
};

void dwdw_colptrs_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 2> dwdw_rowvals_Okuonghae_ChaosSolitonsFractals2020_ = {
    1, 2
};

void dwdw_rowvals_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
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
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void dwdw_Okuonghae_ChaosSolitonsFractals2020(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dlam_dtotal_pop = -transmission_rate_effective*(1 - delta)*(1 - eps)*(alpha*asymptomatic + symptomatic)/std::pow(-detected + total_pop, 2);  // dwdw[0]
    dflux_susceptible_to_exposed_dlam = 1.0*susceptible;  // dwdw[1]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
