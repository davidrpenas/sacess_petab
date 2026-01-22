#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Perelson_Science1996 {

static constexpr std::array<sunindextype, 4> dwdp_colptrs_Perelson_Science1996_ = {
    0, 2, 4, 4
};

void dwdp_colptrs_Perelson_Science1996(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Perelson_Science1996_));
}
} // namespace model_Perelson_Science1996
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Perelson_Science1996 {

static constexpr std::array<sunindextype, 4> dwdp_rowvals_Perelson_Science1996_ = {
    2, 3, 1, 4
};

void dwdp_rowvals_Perelson_Science1996(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Perelson_Science1996_));
}
} // namespace model_Perelson_Science1996
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
namespace model_Perelson_Science1996 {

void dwdp_Perelson_Science1996(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v3_dc = Vin;  // dwdp[0]
    dflux_v4_dc = Vni;  // dwdp[1]
    dflux_v2_ddelta = Tstar;  // dwdp[2]
    dflux_v5_ddelta = NN*Tstar;  // dwdp[3]
}

} // namespace model_Perelson_Science1996
} // namespace amici
