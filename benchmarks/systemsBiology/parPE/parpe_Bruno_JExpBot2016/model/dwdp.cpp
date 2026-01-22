#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 95> dwdp_colptrs_Bruno_JExpBot2016_ = {
    0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12
};

void dwdp_colptrs_Bruno_JExpBot2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 12> dwdp_rowvals_Bruno_JExpBot2016_ = {
    5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4
};

void dwdp_rowvals_Bruno_JExpBot2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "dwdp.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dwdp_Bruno_JExpBot2016(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v6_ReactionName_dk5 = 1.0*k5_multiplier*zea;  // dwdp[0]
    dflux_v1_ReactionName_dkb1 = 1.0*bcar*kb1_multiplier;  // dwdp[1]
    dflux_v2_ReactionName_dkb2 = 1.0*b10*kb2_multiplier;  // dwdp[2]
    dflux_v3_ReactionName_dkc1 = 1.0*bcry*kc1_multiplier;  // dwdp[3]
    dflux_v4_ReactionName_dkc2 = 1.0*bcry*kc2_multiplier;  // dwdp[4]
    dflux_v5_ReactionName_dkc4 = 1.0*kc4_multiplier*ohb10;  // dwdp[5]
    dflux_v6_ReactionName_dk5_multiplier = 1.0*k5*zea;  // dwdp[6]
    dflux_v1_ReactionName_dkb1_multiplier = 1.0*bcar*kb1;  // dwdp[7]
    dflux_v2_ReactionName_dkb2_multiplier = 1.0*b10*kb2;  // dwdp[8]
    dflux_v3_ReactionName_dkc1_multiplier = 1.0*bcry*kc1;  // dwdp[9]
    dflux_v4_ReactionName_dkc2_multiplier = 1.0*bcry*kc2;  // dwdp[10]
    dflux_v5_ReactionName_dkc4_multiplier = 1.0*kc4*ohb10;  // dwdp[11]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
