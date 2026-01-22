#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

static constexpr std::array<sunindextype, 29> dwdp_colptrs_Smith_BMCSystBiol2013_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17
};

void dwdp_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Smith_BMCSystBiol2013_));
}
} // namespace model_Smith_BMCSystBiol2013
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

static constexpr std::array<sunindextype, 17> dwdp_rowvals_Smith_BMCSystBiol2013_ = {
    29, 30, 53, 54, 61, 68, 69, 71, 73, 74, 75, 83, 84, 85, 86, 87, 88
};

void dwdp_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Smith_BMCSystBiol2013_));
}
} // namespace model_Smith_BMCSystBiol2013
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
#include "dwdp.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void dwdp_Smith_BMCSystBiol2013(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_R1f_dk1 = InR*Ins;  // dwdp[0]
    dflux_R1r_dkminus1 = Ins_InR;  // dwdp[1]
    dflux_R30f_dk30f = PTP1B*ROS;  // dwdp[2]
    dflux_R30r_dk30r = GSH*PTP1B_ox;  // dwdp[3]
    dflux_R35f_dk35f = NOX;  // dwdp[4]
    dflux_R12f_dk7 = IRS1*cyto_vol*(Ins_2_InR_P + Ins_InR_P)/IRp;  // dwdp[5]
    dflux_R12r_dkminus7a = IRS1_TyrP*PTP1B*cyto_vol;  // dwdp[6]
    dflux_R12_a_r_dkminus7b = IRS1_PolySerP*PP2A*cyto_vol;  // dwdp[7]
    dflux_R12_b_r_dkminus7b = IRS1_TyrP_PolySerP*PP2A*cyto_vol;  // dwdp[8]
    dflux_R13f_dk8 = IRS1_TyrP*PI3K*cyto_vol;  // dwdp[9]
    dflux_R13r_dkminus8 = IRS1_TyrP_PI3K*cyto_vol;  // dwdp[10]
    dflux_R42f_dk42f = JNK*ROS*alpha_ox*indicator_jnk;  // dwdp[11]
    dflux_R42r_dk42r = DUSP*JNK_P*indicator_jnk;  // dwdp[12]
    dflux_R43f_dk43f = IKK*ROS*indicator_jnk;  // dwdp[13]
    dflux_R43r_dk43r = DUSP*IKK_P*indicator_jnk;  // dwdp[14]
    dflux_R32f_dk32f = DUSP*ROS*indicator_jnk;  // dwdp[15]
    dflux_R32r_dk32r = DUSP_ox*GSH*indicator_jnk;  // dwdp[16]
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
