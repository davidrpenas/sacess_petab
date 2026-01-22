#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Perelson_Science1996 {

static constexpr std::array<sunindextype, 6> dxdotdw_colptrs_Perelson_Science1996_ = {
    0, 1, 2, 4, 6, 8
};

void dxdotdw_colptrs_Perelson_Science1996(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Perelson_Science1996_));
}
} // namespace model_Perelson_Science1996
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Perelson_Science1996 {

static constexpr std::array<sunindextype, 8> dxdotdw_rowvals_Perelson_Science1996_ = {
    0, 0, 1, 2, 1, 3, 1, 3
};

void dxdotdw_rowvals_Perelson_Science1996(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Perelson_Science1996_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Perelson_Science1996 {

void dxdotdw_Perelson_Science1996(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddTstardt_dflux_v1 = 1;  // dxdotdw[0]
    ddTstardt_dflux_v2 = -1;  // dxdotdw[1]
    ddVdt_dflux_v3 = -1;  // dxdotdw[2]
    ddVindt_dflux_v3 = -1;  // dxdotdw[3]
    ddVdt_dflux_v4 = -1;  // dxdotdw[4]
    ddVnidt_dflux_v4 = -1;  // dxdotdw[5]
    ddVdt_dflux_v5 = 1;  // dxdotdw[6]
    ddVnidt_dflux_v5 = 1;  // dxdotdw[7]
}

} // namespace model_Perelson_Science1996
} // namespace amici
