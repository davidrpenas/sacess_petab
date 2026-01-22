#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 27> dxdotdw_colptrs_Weber_BMC2015_ = {
    0, 0, 0, 0, 0, 0, 2, 4, 6, 7, 8, 9, 10, 12, 14, 15, 16, 17, 18, 20, 22, 24, 25, 26, 27, 28, 29
};

void dxdotdw_colptrs_Weber_BMC2015(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Weber_BMC2015_));
}
} // namespace model_Weber_BMC2015
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 29> dxdotdw_rowvals_Weber_BMC2015_ = {
    0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 2, 3, 2, 3, 2, 2, 2, 3, 4, 6, 4, 5, 5, 6, 4, 4, 4, 5, 6
};

void dxdotdw_rowvals_Weber_BMC2015(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Weber_BMC2015_));
}
} // namespace model_Weber_BMC2015
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
#include "dxdotdw.h"

namespace amici {
namespace model_Weber_BMC2015 {

void dxdotdw_Weber_BMC2015(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddPKDdt_dflux_v1_v_0 = -1;  // dxdotdw[0]
    ddPKDDAGadt_dflux_v1_v_0 = 1;  // dxdotdw[1]
    ddPKDdt_dflux_v2_v_1 = -1;  // dxdotdw[2]
    ddPKDDAGadt_dflux_v2_v_1 = 1;  // dxdotdw[3]
    ddPKDdt_dflux_v3_v_2 = 1;  // dxdotdw[4]
    ddPKDDAGadt_dflux_v3_v_2 = -1;  // dxdotdw[5]
    ddPKDdt_dflux_v4_v_3 = 1;  // dxdotdw[6]
    ddPKDdt_dflux_v5_v_4 = 1;  // dxdotdw[7]
    ddPKDdt_dflux_v6_v_5 = -1;  // dxdotdw[8]
    ddPKDDAGadt_dflux_v7_v_6 = -1;  // dxdotdw[9]
    ddPI4K3Bdt_dflux_v8_v_7 = 1;  // dxdotdw[10]
    ddPI4K3Badt_dflux_v8_v_7 = -1;  // dxdotdw[11]
    ddPI4K3Bdt_dflux_v9_v_8 = -1;  // dxdotdw[12]
    ddPI4K3Badt_dflux_v9_v_8 = 1;  // dxdotdw[13]
    ddPI4K3Bdt_dflux_v10_v_9 = 1;  // dxdotdw[14]
    ddPI4K3Bdt_dflux_v11_v_10 = 1;  // dxdotdw[15]
    ddPI4K3Bdt_dflux_v12_v_11 = -1;  // dxdotdw[16]
    ddPI4K3Badt_dflux_v13_v_12 = -1;  // dxdotdw[17]
    ddCERTERadt_dflux_v14_v_13 = -1;  // dxdotdw[18]
    ddCERTTGNadt_dflux_v14_v_13 = 1;  // dxdotdw[19]
    ddCERTERadt_dflux_v15_v_14 = 1;  // dxdotdw[20]
    ddCERTdt_dflux_v15_v_14 = -1;  // dxdotdw[21]
    ddCERTdt_dflux_v16_v_15 = 1;  // dxdotdw[22]
    ddCERTTGNadt_dflux_v16_v_15 = -1;  // dxdotdw[23]
    ddCERTERadt_dflux_v17_v_16 = 1;  // dxdotdw[24]
    ddCERTERadt_dflux_v18_v_17 = 1;  // dxdotdw[25]
    ddCERTERadt_dflux_v19_v_18 = -1;  // dxdotdw[26]
    ddCERTdt_dflux_v20_v_19 = -1;  // dxdotdw[27]
    ddCERTTGNadt_dflux_v21_v_20 = -1;  // dxdotdw[28]
}

} // namespace model_Weber_BMC2015
} // namespace amici
