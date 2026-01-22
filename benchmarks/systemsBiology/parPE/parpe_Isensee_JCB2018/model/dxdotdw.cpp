#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 52> dxdotdw_colptrs_Isensee_JCB2018_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 22, 23, 24, 27, 30, 33, 36, 38, 40, 42, 44, 46, 48, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 74, 76, 78, 80, 82, 84
};

void dxdotdw_colptrs_Isensee_JCB2018(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 84> dxdotdw_rowvals_Isensee_JCB2018_ = {
    11, 12, 13, 14, 22, 0, 2, 0, 2, 1, 3, 1, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 4, 4, 5, 7, 21, 5, 7, 21, 5, 8, 20, 5, 8, 20, 20, 21, 7, 8, 7, 8, 8, 9, 8, 9, 8, 9, 5, 9, 10, 10, 20, 10, 20, 8, 15, 8, 15, 8, 16, 8, 16, 8, 17, 8, 17, 8, 18, 8, 18, 5, 18, 19, 19, 20, 5, 6, 5, 6, 23, 24, 23, 24
};

void dxdotdw_rowvals_Isensee_JCB2018(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
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
namespace model_Isensee_JCB2018 {

void dxdotdw_Isensee_JCB2018(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSp8_Br_cAMPSdt_dflux_v1_v_0 = 1;  // dxdotdw[0]
    ddRp8_Br_cAMPSdt_dflux_v2_v_1 = 1;  // dxdotdw[1]
    ddRp8_pCPT_cAMPSdt_dflux_v3_v_2 = 1;  // dxdotdw[2]
    ddRp_cAMPSdt_dflux_v4_v_3 = 1;  // dxdotdw[3]
    ddIBMXdt_dflux_v5_v_4 = 1;  // dxdotdw[4]
    ddACdt_dflux_v6_v_5 = -1;  // dxdotdw[5]
    ddAC_Fskdt_dflux_v6_v_5 = 1;  // dxdotdw[6]
    ddACdt_dflux_v7_v_6 = 1;  // dxdotdw[7]
    ddAC_Fskdt_dflux_v7_v_6 = -1;  // dxdotdw[8]
    ddpACdt_dflux_v8_v_7 = -1;  // dxdotdw[9]
    ddpAC_Fskdt_dflux_v8_v_7 = 1;  // dxdotdw[10]
    ddpACdt_dflux_v9_v_8 = 1;  // dxdotdw[11]
    ddpAC_Fskdt_dflux_v9_v_8 = -1;  // dxdotdw[12]
    ddACdt_dflux_v10_v_9 = -1;  // dxdotdw[13]
    ddpACdt_dflux_v10_v_9 = 1;  // dxdotdw[14]
    ddAC_Fskdt_dflux_v11_v_10 = -1;  // dxdotdw[15]
    ddpAC_Fskdt_dflux_v11_v_10 = 1;  // dxdotdw[16]
    ddACdt_dflux_v12_v_11 = 1;  // dxdotdw[17]
    ddpACdt_dflux_v12_v_11 = -1;  // dxdotdw[18]
    ddAC_Fskdt_dflux_v13_v_12 = 1;  // dxdotdw[19]
    ddpAC_Fskdt_dflux_v13_v_12 = -1;  // dxdotdw[20]
    ddcAMPdt_dflux_v14_v_13 = 1;  // dxdotdw[21]
    ddcAMPdt_dflux_v15_v_14 = 1;  // dxdotdw[22]
    ddcAMPdt_dflux_v16_v_15 = -1;  // dxdotdw[23]
    ddCsubdt_dflux_v17_v_16 = 1;  // dxdotdw[24]
    ddRII_C_2dt_dflux_v17_v_16 = -1;  // dxdotdw[25]
    ddRII_2dt_dflux_v17_v_16 = 1;  // dxdotdw[26]
    ddCsubdt_dflux_v18_v_17 = -1;  // dxdotdw[27]
    ddRII_C_2dt_dflux_v18_v_17 = 1;  // dxdotdw[28]
    ddRII_2dt_dflux_v18_v_17 = -1;  // dxdotdw[29]
    ddCsubdt_dflux_v19_v_18 = 1;  // dxdotdw[30]
    ddRIIp_C_2dt_dflux_v19_v_18 = -1;  // dxdotdw[31]
    ddRIIp_2dt_dflux_v19_v_18 = 1;  // dxdotdw[32]
    ddCsubdt_dflux_v20_v_19 = -1;  // dxdotdw[33]
    ddRIIp_C_2dt_dflux_v20_v_19 = 1;  // dxdotdw[34]
    ddRIIp_2dt_dflux_v20_v_19 = -1;  // dxdotdw[35]
    ddRIIp_2dt_dflux_v21_v_20 = -1;  // dxdotdw[36]
    ddRII_2dt_dflux_v21_v_20 = 1;  // dxdotdw[37]
    ddRII_C_2dt_dflux_v22_v_21 = -1;  // dxdotdw[38]
    ddRIIp_C_2dt_dflux_v22_v_21 = 1;  // dxdotdw[39]
    ddRII_C_2dt_dflux_v23_v_22 = 1;  // dxdotdw[40]
    ddRIIp_C_2dt_dflux_v23_v_22 = -1;  // dxdotdw[41]
    ddRIIp_C_2dt_dflux_v24_v_23 = -1;  // dxdotdw[42]
    ddRIIp_cAMP_C_2dt_dflux_v24_v_23 = 1;  // dxdotdw[43]
    ddRIIp_C_2dt_dflux_v25_v_24 = 1;  // dxdotdw[44]
    ddRIIp_cAMP_C_2dt_dflux_v25_v_24 = -1;  // dxdotdw[45]
    ddRIIp_C_2dt_dflux_v26_v_25 = 1;  // dxdotdw[46]
    ddRIIp_cAMP_C_2dt_dflux_v26_v_25 = -1;  // dxdotdw[47]
    ddCsubdt_dflux_v27_v_26 = 1;  // dxdotdw[48]
    ddRIIp_cAMP_C_2dt_dflux_v27_v_26 = -1;  // dxdotdw[49]
    ddRIIp_cAMP_2dt_dflux_v27_v_26 = 1;  // dxdotdw[50]
    ddRIIp_cAMP_2dt_dflux_v28_v_27 = -1;  // dxdotdw[51]
    ddRIIp_2dt_dflux_v28_v_27 = 1;  // dxdotdw[52]
    ddRIIp_cAMP_2dt_dflux_v29_v_28 = -1;  // dxdotdw[53]
    ddRIIp_2dt_dflux_v29_v_28 = 1;  // dxdotdw[54]
    ddRIIp_C_2dt_dflux_v30_v_29 = -1;  // dxdotdw[55]
    ddRIIp_Rp8_Br_cAMPS_C_2dt_dflux_v30_v_29 = 1;  // dxdotdw[56]
    ddRIIp_C_2dt_dflux_v31_v_30 = 1;  // dxdotdw[57]
    ddRIIp_Rp8_Br_cAMPS_C_2dt_dflux_v31_v_30 = -1;  // dxdotdw[58]
    ddRIIp_C_2dt_dflux_v32_v_31 = -1;  // dxdotdw[59]
    ddRIIp_Rp8_pCPT_cAMPS_C_2dt_dflux_v32_v_31 = 1;  // dxdotdw[60]
    ddRIIp_C_2dt_dflux_v33_v_32 = 1;  // dxdotdw[61]
    ddRIIp_Rp8_pCPT_cAMPS_C_2dt_dflux_v33_v_32 = -1;  // dxdotdw[62]
    ddRIIp_C_2dt_dflux_v34_v_33 = -1;  // dxdotdw[63]
    ddRIIp_Rp_cAMPS_C_2dt_dflux_v34_v_33 = 1;  // dxdotdw[64]
    ddRIIp_C_2dt_dflux_v35_v_34 = 1;  // dxdotdw[65]
    ddRIIp_Rp_cAMPS_C_2dt_dflux_v35_v_34 = -1;  // dxdotdw[66]
    ddRIIp_C_2dt_dflux_v36_v_35 = -1;  // dxdotdw[67]
    ddRIIp_Sp8_Br_cAMPS_C_2dt_dflux_v36_v_35 = 1;  // dxdotdw[68]
    ddRIIp_C_2dt_dflux_v37_v_36 = 1;  // dxdotdw[69]
    ddRIIp_Sp8_Br_cAMPS_C_2dt_dflux_v37_v_36 = -1;  // dxdotdw[70]
    ddCsubdt_dflux_v38_v_37 = 1;  // dxdotdw[71]
    ddRIIp_Sp8_Br_cAMPS_C_2dt_dflux_v38_v_37 = -1;  // dxdotdw[72]
    ddRIIp_Sp8_Br_cAMPS_2dt_dflux_v38_v_37 = 1;  // dxdotdw[73]
    ddRIIp_Sp8_Br_cAMPS_2dt_dflux_v39_v_38 = -1;  // dxdotdw[74]
    ddRIIp_2dt_dflux_v39_v_38 = 1;  // dxdotdw[75]
    ddCsubdt_dflux_v40_v_39 = -1;  // dxdotdw[76]
    ddCsub_H89dt_dflux_v40_v_39 = 1;  // dxdotdw[77]
    ddCsubdt_dflux_v41_v_40 = 1;  // dxdotdw[78]
    ddCsub_H89dt_dflux_v41_v_40 = -1;  // dxdotdw[79]
    ddPDEdt_dflux_v42_v_41 = -1;  // dxdotdw[80]
    ddpPDEdt_dflux_v42_v_41 = 1;  // dxdotdw[81]
    ddPDEdt_dflux_v43_v_42 = 1;  // dxdotdw[82]
    ddpPDEdt_dflux_v43_v_42 = -1;  // dxdotdw[83]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
