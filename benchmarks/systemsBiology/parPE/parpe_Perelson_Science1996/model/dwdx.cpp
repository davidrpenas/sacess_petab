#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Perelson_Science1996 {

static constexpr std::array<sunindextype, 5> dwdx_colptrs_Perelson_Science1996_ = {
    0, 2, 2, 4, 5
};

void dwdx_colptrs_Perelson_Science1996(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Perelson_Science1996_));
}
} // namespace model_Perelson_Science1996
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Perelson_Science1996 {

static constexpr std::array<sunindextype, 5> dwdx_rowvals_Perelson_Science1996_ = {
    1, 4, 0, 2, 3
};

void dwdx_rowvals_Perelson_Science1996(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Perelson_Science1996_));
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
#include "dwdx.h"

namespace amici {
namespace model_Perelson_Science1996 {

void dwdx_Perelson_Science1996(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v2_dTstar = delta;  // dwdx[0]
        dflux_v5_dTstar = NN*delta;  // dwdx[1]
        dflux_v1_dVin = K0*T0;  // dwdx[2]
        dflux_v3_dVin = c;  // dwdx[3]
        dflux_v4_dVni = c;  // dwdx[4]
    }
}

} // namespace model_Perelson_Science1996
} // namespace amici
