#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<sunindextype, 13> dxdotdw_colptrs_Fujita_SciSignal2010_ = {
    0, 0, 2, 5, 8, 9, 12, 15, 17, 19, 21, 22, 23
};

void dxdotdw_colptrs_Fujita_SciSignal2010(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Fujita_SciSignal2010_));
}
} // namespace model_Fujita_SciSignal2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<sunindextype, 23> dxdotdw_rowvals_Fujita_SciSignal2010_ = {
    0, 8, 1, 2, 3, 1, 2, 4, 1, 4, 5, 6, 4, 6, 7, 3, 4, 5, 7, 1, 8, 0, 0
};

void dxdotdw_rowvals_Fujita_SciSignal2010(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Fujita_SciSignal2010_));
}
} // namespace model_Fujita_SciSignal2010
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
namespace model_Fujita_SciSignal2010 {

void dxdotdw_Fujita_SciSignal2010(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddEGFRdt_dflux_v1_reaction_1 = -1;  // dxdotdw[0]
    ddEGF_EGFRdt_dflux_v1_reaction_1 = 1;  // dxdotdw[1]
    ddpEGFRdt_dflux_v2_reaction_2 = -1;  // dxdotdw[2]
    ddpEGFR_Aktdt_dflux_v2_reaction_2 = 1;  // dxdotdw[3]
    ddAktdt_dflux_v2_reaction_2 = -1;  // dxdotdw[4]
    ddpEGFRdt_dflux_v3_reaction_3 = 1;  // dxdotdw[5]
    ddpEGFR_Aktdt_dflux_v3_reaction_3 = -1;  // dxdotdw[6]
    ddpAktdt_dflux_v3_reaction_3 = 1;  // dxdotdw[7]
    ddpEGFRdt_dflux_v4_reaction_4 = -1;  // dxdotdw[8]
    ddpAktdt_dflux_v5_reaction_5 = -1;  // dxdotdw[9]
    ddS6dt_dflux_v5_reaction_5 = -1;  // dxdotdw[10]
    ddpAkt_S6dt_dflux_v5_reaction_5 = 1;  // dxdotdw[11]
    ddpAktdt_dflux_v6_reaction_6 = 1;  // dxdotdw[12]
    ddpAkt_S6dt_dflux_v6_reaction_6 = -1;  // dxdotdw[13]
    ddpS6dt_dflux_v6_reaction_6 = 1;  // dxdotdw[14]
    ddAktdt_dflux_v7_reaction_7 = 1;  // dxdotdw[15]
    ddpAktdt_dflux_v7_reaction_7 = -1;  // dxdotdw[16]
    ddS6dt_dflux_v8_reaction_8 = 1;  // dxdotdw[17]
    ddpS6dt_dflux_v8_reaction_8 = -1;  // dxdotdw[18]
    ddpEGFRdt_dflux_v9_reaction_9 = 1;  // dxdotdw[19]
    ddEGF_EGFRdt_dflux_v9_reaction_9 = -1;  // dxdotdw[20]
    ddEGFRdt_dflux_v10_reaction_10 = -1;  // dxdotdw[21]
    ddEGFRdt_dflux_v11_reaction_11 = 1;  // dxdotdw[22]
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
