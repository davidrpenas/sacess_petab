#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<sunindextype, 17> dxdotdw_colptrs_Raia_CancerResearch2011_ = {
    0, 2, 4, 6, 8, 10, 12, 14, 15, 17, 19, 21, 23, 24, 25, 26, 27
};

void dxdotdw_colptrs_Raia_CancerResearch2011(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Raia_CancerResearch2011_));
}
} // namespace model_Raia_CancerResearch2011
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<sunindextype, 27> dxdotdw_rowvals_Raia_CancerResearch2011_ = {
    0, 2, 0, 1, 0, 1, 2, 3, 5, 6, 5, 6, 3, 4, 4, 5, 6, 7, 8, 7, 8, 10, 11, 9, 12, 12, 13
};

void dxdotdw_rowvals_Raia_CancerResearch2011(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Raia_CancerResearch2011_));
}
} // namespace model_Raia_CancerResearch2011
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
namespace model_Raia_CancerResearch2011 {

void dxdotdw_Raia_CancerResearch2011(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddRecdt_dflux_v1_v_0 = -0.01;  // dxdotdw[0]
    ddIL13_Recdt_dflux_v1_v_0 = 0.01;  // dxdotdw[1]
    ddRecdt_dflux_v2_v_1 = -0.01;  // dxdotdw[2]
    ddRec_idt_dflux_v2_v_1 = 0.01;  // dxdotdw[3]
    ddRecdt_dflux_v3_v_2 = 0.01;  // dxdotdw[4]
    ddRec_idt_dflux_v3_v_2 = -0.01;  // dxdotdw[5]
    ddIL13_Recdt_dflux_v4_v_3 = -0.01;  // dxdotdw[6]
    ddp_IL13_Recdt_dflux_v4_v_3 = 0.01;  // dxdotdw[7]
    ddJAK2dt_dflux_v5_v_4 = -0.01;  // dxdotdw[8]
    ddpJAK2dt_dflux_v5_v_4 = 0.01;  // dxdotdw[9]
    ddJAK2dt_dflux_v6_v_5 = -0.01;  // dxdotdw[10]
    ddpJAK2dt_dflux_v6_v_5 = 0.01;  // dxdotdw[11]
    ddp_IL13_Recdt_dflux_v7_v_6 = -0.01;  // dxdotdw[12]
    ddp_IL13_Rec_idt_dflux_v7_v_6 = 0.01;  // dxdotdw[13]
    ddp_IL13_Rec_idt_dflux_v8_v_7 = -0.01;  // dxdotdw[14]
    ddJAK2dt_dflux_v9_v_8 = 0.01;  // dxdotdw[15]
    ddpJAK2dt_dflux_v9_v_8 = -0.01;  // dxdotdw[16]
    ddSTAT5dt_dflux_v10_v_9 = -0.01;  // dxdotdw[17]
    ddpSTAT5dt_dflux_v10_v_9 = 0.01;  // dxdotdw[18]
    ddSTAT5dt_dflux_v11_v_10 = 0.01;  // dxdotdw[19]
    ddpSTAT5dt_dflux_v11_v_10 = -0.01;  // dxdotdw[20]
    ddDecoyRdt_dflux_v12_v_11 = -0.01;  // dxdotdw[21]
    ddIL13_DecoyRdt_dflux_v12_v_11 = 0.01;  // dxdotdw[22]
    ddSOCS3mRNAdt_dflux_v13_v_12 = 0.10000000000000001;  // dxdotdw[23]
    ddSOCS3dt_dflux_v14_v_13 = 0.01;  // dxdotdw[24]
    ddSOCS3dt_dflux_v15_v_14 = -0.01;  // dxdotdw[25]
    ddCD274mRNAdt_dflux_v16_v_15 = 0.10000000000000001;  // dxdotdw[26]
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
