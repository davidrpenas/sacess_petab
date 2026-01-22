#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 14> dxdotdw_colptrs_Sneyd_PNAS2002_ = {
    0, 0, 0, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20
};

void dxdotdw_colptrs_Sneyd_PNAS2002(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Sneyd_PNAS2002_));
}
} // namespace model_Sneyd_PNAS2002
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 20> dxdotdw_rowvals_Sneyd_PNAS2002_ = {
    0, 1, 0, 1, 1, 2, 1, 2, 0, 4, 0, 4, 4, 5, 4, 5, 0, 3, 0, 3
};

void dxdotdw_rowvals_Sneyd_PNAS2002(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Sneyd_PNAS2002_));
}
} // namespace model_Sneyd_PNAS2002
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
namespace model_Sneyd_PNAS2002 {

void dxdotdw_Sneyd_PNAS2002(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddIPR_Odt_dflux_v_0 = -1;  // dxdotdw[0]
    ddIPR_Rdt_dflux_v_0 = 1;  // dxdotdw[1]
    ddIPR_Odt_dflux_v_1 = 1;  // dxdotdw[2]
    ddIPR_Rdt_dflux_v_1 = -1;  // dxdotdw[3]
    ddIPR_Rdt_dflux_v_2 = -1;  // dxdotdw[4]
    ddIPR_I1dt_dflux_v_2 = 1;  // dxdotdw[5]
    ddIPR_Rdt_dflux_v_3 = 1;  // dxdotdw[6]
    ddIPR_I1dt_dflux_v_3 = -1;  // dxdotdw[7]
    ddIPR_Odt_dflux_v_4 = -1;  // dxdotdw[8]
    ddIPR_Adt_dflux_v_4 = 1;  // dxdotdw[9]
    ddIPR_Odt_dflux_v_5 = 1;  // dxdotdw[10]
    ddIPR_Adt_dflux_v_5 = -1;  // dxdotdw[11]
    ddIPR_Adt_dflux_v_6 = -1;  // dxdotdw[12]
    ddIPR_I2dt_dflux_v_6 = 1;  // dxdotdw[13]
    ddIPR_Adt_dflux_v_7 = 1;  // dxdotdw[14]
    ddIPR_I2dt_dflux_v_7 = -1;  // dxdotdw[15]
    ddIPR_Odt_dflux_v_8 = -1;  // dxdotdw[16]
    ddIPR_Sdt_dflux_v_8 = 1;  // dxdotdw[17]
    ddIPR_Odt_dflux_v_9 = 1;  // dxdotdw[18]
    ddIPR_Sdt_dflux_v_9 = -1;  // dxdotdw[19]
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
