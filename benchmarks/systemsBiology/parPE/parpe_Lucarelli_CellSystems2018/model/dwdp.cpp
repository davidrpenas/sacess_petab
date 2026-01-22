#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<sunindextype, 85> dwdp_colptrs_Lucarelli_CellSystems2018_ = {
    0, 1, 1, 1, 1, 3, 17, 19, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 56, 57, 58, 59, 60, 62, 63, 64, 65, 66, 67, 69, 70, 71, 72, 74, 75, 76, 77, 79, 80, 81, 82, 83, 85, 86, 87, 88, 89, 91, 91, 92, 93, 94, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95
};

void dwdp_colptrs_Lucarelli_CellSystems2018(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Lucarelli_CellSystems2018_));
}
} // namespace model_Lucarelli_CellSystems2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<sunindextype, 95> dwdp_rowvals_Lucarelli_CellSystems2018_ = {
    0, 12, 15, 5, 7, 11, 14, 17, 18, 20, 22, 23, 25, 27, 29, 31, 32, 10, 13, 33, 33, 33, 33, 33, 34, 35, 35, 35, 35, 35, 35, 36, 37, 37, 37, 37, 38, 39, 39, 39, 39, 39, 40, 41, 41, 41, 41, 41, 41, 42, 43, 43, 43, 43, 43, 44, 45, 45, 45, 45, 45, 46, 47, 47, 47, 47, 47, 47, 48, 49, 49, 49, 49, 50, 51, 51, 51, 51, 52, 53, 53, 53, 53, 53, 54, 55, 55, 55, 55, 55, 56, 21, 30, 26, 1
};

void dwdp_rowvals_Lucarelli_CellSystems2018(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Lucarelli_CellSystems2018_));
}
} // namespace model_Lucarelli_CellSystems2018
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
namespace model_Lucarelli_CellSystems2018 {

void dwdp_Lucarelli_CellSystems2018(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v1_v_0_dRec_act = 13.17*Rec*TGFb;  // dwdp[0]
    dflux_v13_v_12_dS_dephos = 13.17*pS2;  // dwdp[1]
    dflux_v16_v_15_dS_dephos = 13.17*pS3;  // dwdp[2]
    dflux_v6_v_5_dS_dephosphos = 39.509999999999998*ppS2_ppS2_ppS2;  // dwdp[3]
    dflux_v8_v_7_dS_dephosphos = 39.509999999999998*ppS3_ppS3_ppS3;  // dwdp[4]
    dflux_v12_v_11_dS_dephosphos = 13.17*ppS2;  // dwdp[5]
    dflux_v15_v_14_dS_dephosphos = 13.17*ppS3;  // dwdp[6]
    dflux_v18_v_17_dS_dephosphos = 26.34*ppS2_ppS2_ppS3;  // dwdp[7]
    dflux_v19_v_18_dS_dephosphos = 13.17*ppS2_ppS2_ppS3;  // dwdp[8]
    dflux_v21_v_20_dS_dephosphos = 26.34*ppS2_ppS2_S4;  // dwdp[9]
    dflux_v23_v_22_dS_dephosphos = 13.17*ppS2_ppS3_ppS3;  // dwdp[10]
    dflux_v24_v_23_dS_dephosphos = 26.34*ppS2_ppS3_ppS3;  // dwdp[11]
    dflux_v26_v_25_dS_dephosphos = 26.34*ppS3_ppS3_S4;  // dwdp[12]
    dflux_v28_v_27_dS_dephosphos = 13.17*ppS2_S4_S4;  // dwdp[13]
    dflux_v30_v_29_dS_dephosphos = 13.17*ppS3_S4_S4;  // dwdp[14]
    dflux_v32_v_31_dS_dephosphos = 13.17*ppS2_ppS3_S4;  // dwdp[15]
    dflux_v33_v_32_dS_dephosphos = 13.17*ppS2_ppS3_S4;  // dwdp[16]
    dflux_v11_v_10_dS_phos = 13.17*S2*TGFb_pRec;  // dwdp[17]
    dflux_v14_v_13_dS_phos = 13.17*S3*TGFb_pRec;  // dwdp[18]
    dflux_v34_v_33_dgeneA_act1 = 13.17*ppS2_ppS3_S4/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[19]
    dflux_v34_v_33_dgeneA_act2 = 13.17*ppS2_S4_S4/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[20]
    dflux_v34_v_33_dgeneA_inh2 = -ppS2_S4_S4*(13.17*geneA_act1*ppS2_ppS3_S4 + 13.17*geneA_act2*ppS2_S4_S4 + 13.17*geneA_act3*ppS2_ppS3_ppS3 + 13.17*geneA_turn)/std::pow(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[21]
    dflux_v34_v_33_dgeneA_inh3 = -ppS2_ppS3_ppS3*(13.17*geneA_act1*ppS2_ppS3_S4 + 13.17*geneA_act2*ppS2_S4_S4 + 13.17*geneA_act3*ppS2_ppS3_ppS3 + 13.17*geneA_turn)/std::pow(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[22]
    dflux_v34_v_33_dgeneA_turn = 13.17/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[23]
    dflux_v35_v_34_dgeneA_turn = 13.17*geneA;  // dwdp[24]
    dflux_v36_v_35_dgeneB_act1 = 13.17*ppS2_ppS3_S4/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[25]
    dflux_v36_v_35_dgeneB_act2 = 13.17*ppS2_S4_S4/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[26]
    dflux_v36_v_35_dgeneB_act3 = 13.17*ppS2_ppS3_ppS3/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[27]
    dflux_v36_v_35_dgeneB_inh2 = -ppS2_S4_S4*(13.17*geneB_act1*ppS2_ppS3_S4 + 13.17*geneB_act2*ppS2_S4_S4 + 13.17*geneB_act3*ppS2_ppS3_ppS3 + 13.17*geneB_turn)/std::pow(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[28]
    dflux_v36_v_35_dgeneB_inh3 = -ppS2_ppS3_ppS3*(13.17*geneB_act1*ppS2_ppS3_S4 + 13.17*geneB_act2*ppS2_S4_S4 + 13.17*geneB_act3*ppS2_ppS3_ppS3 + 13.17*geneB_turn)/std::pow(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[29]
    dflux_v36_v_35_dgeneB_turn = 13.17/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[30]
    dflux_v37_v_36_dgeneB_turn = 13.17*geneB;  // dwdp[31]
    dflux_v38_v_37_dgeneC_act1 = 13.17*ppS2_ppS3_S4/(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[32]
    dflux_v38_v_37_dgeneC_inh2 = -ppS2_S4_S4*(13.17*geneC_act1*ppS2_ppS3_S4 + 13.17*geneC_act2*ppS2_S4_S4 + 13.17*geneC_act3*ppS2_ppS3_ppS3 + 13.17*geneC_turn)/std::pow(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[33]
    dflux_v38_v_37_dgeneC_inh3 = -ppS2_ppS3_ppS3*(13.17*geneC_act1*ppS2_ppS3_S4 + 13.17*geneC_act2*ppS2_S4_S4 + 13.17*geneC_act3*ppS2_ppS3_ppS3 + 13.17*geneC_turn)/std::pow(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[34]
    dflux_v38_v_37_dgeneC_turn = 13.17/(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[35]
    dflux_v39_v_38_dgeneC_turn = 13.17*geneC;  // dwdp[36]
    dflux_v40_v_39_dgeneD_act1 = 13.17*ppS2_ppS3_S4/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[37]
    dflux_v40_v_39_dgeneD_act2 = 13.17*ppS2_S4_S4/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[38]
    dflux_v40_v_39_dgeneD_inh2 = -ppS2_S4_S4*(13.17*geneD_act1*ppS2_ppS3_S4 + 13.17*geneD_act2*ppS2_S4_S4 + 13.17*geneD_act3*ppS2_ppS3_ppS3 + 13.17*geneD_turn)/std::pow(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[39]
    dflux_v40_v_39_dgeneD_inh3 = -ppS2_ppS3_ppS3*(13.17*geneD_act1*ppS2_ppS3_S4 + 13.17*geneD_act2*ppS2_S4_S4 + 13.17*geneD_act3*ppS2_ppS3_ppS3 + 13.17*geneD_turn)/std::pow(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[40]
    dflux_v40_v_39_dgeneD_turn = 13.17/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[41]
    dflux_v41_v_40_dgeneD_turn = 13.17*geneD;  // dwdp[42]
    dflux_v42_v_41_dgeneE_act2 = 13.17*ppS2_S4_S4/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[43]
    dflux_v42_v_41_dgeneE_act3 = 13.17*ppS2_ppS3_ppS3/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[44]
    dflux_v42_v_41_dgeneE_inh1 = -ppS2_ppS3_S4*(13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/std::pow(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[45]
    dflux_v42_v_41_dgeneE_inh2 = -ppS2_S4_S4*(13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/std::pow(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[46]
    dflux_v42_v_41_dgeneE_inh3 = -ppS2_ppS3_ppS3*(13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/std::pow(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[47]
    dflux_v42_v_41_dgeneE_turn = 13.17/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[48]
    dflux_v43_v_42_dgeneE_turn = 13.17*geneE;  // dwdp[49]
    dflux_v44_v_43_dgeneF_act1 = 13.17*ppS2_ppS3_S4/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[50]
    dflux_v44_v_43_dgeneF_act2 = 13.17*ppS2_S4_S4/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[51]
    dflux_v44_v_43_dgeneF_act3 = 13.17*ppS2_ppS3_ppS3/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[52]
    dflux_v44_v_43_dgeneF_inh1 = -ppS2_ppS3_S4*(13.17*geneF_act1*ppS2_ppS3_S4 + 13.17*geneF_act2*ppS2_S4_S4 + 13.17*geneF_act3*ppS2_ppS3_ppS3 + 13.17*geneF_turn)/std::pow(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[53]
    dflux_v44_v_43_dgeneF_turn = 13.17/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[54]
    dflux_v45_v_44_dgeneF_turn = 13.17*geneF;  // dwdp[55]
    dflux_v46_v_45_dgeneG_act1 = 13.17*ppS2_ppS3_S4/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[56]
    dflux_v46_v_45_dgeneG_act2 = 13.17*ppS2_S4_S4/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[57]
    dflux_v46_v_45_dgeneG_inh2 = -ppS2_S4_S4*(13.17*geneG_act1*ppS2_ppS3_S4 + 13.17*geneG_act2*ppS2_S4_S4 + 13.17*geneG_act3*ppS2_ppS3_ppS3 + 13.17*geneG_turn)/std::pow(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[58]
    dflux_v46_v_45_dgeneG_inh3 = -ppS2_ppS3_ppS3*(13.17*geneG_act1*ppS2_ppS3_S4 + 13.17*geneG_act2*ppS2_S4_S4 + 13.17*geneG_act3*ppS2_ppS3_ppS3 + 13.17*geneG_turn)/std::pow(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[59]
    dflux_v46_v_45_dgeneG_turn = 13.17/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[60]
    dflux_v47_v_46_dgeneG_turn = 13.17*geneG;  // dwdp[61]
    dflux_v48_v_47_dgeneH_act2 = 13.17*ppS2_S4_S4/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[62]
    dflux_v48_v_47_dgeneH_act3 = 13.17*ppS2_ppS3_ppS3/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[63]
    dflux_v48_v_47_dgeneH_inh1 = -ppS2_ppS3_S4*(13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/std::pow(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[64]
    dflux_v48_v_47_dgeneH_inh2 = -ppS2_S4_S4*(13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/std::pow(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[65]
    dflux_v48_v_47_dgeneH_inh3 = -ppS2_ppS3_ppS3*(13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/std::pow(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[66]
    dflux_v48_v_47_dgeneH_turn = 13.17/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[67]
    dflux_v49_v_48_dgeneH_turn = 13.17*geneH;  // dwdp[68]
    dflux_v50_v_49_dgeneI_act1 = 13.17*ppS2_ppS3_S4/(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[69]
    dflux_v50_v_49_dgeneI_inh2 = -ppS2_S4_S4*(13.17*geneI_act1*ppS2_ppS3_S4 + 13.17*geneI_act2*ppS2_S4_S4 + 13.17*geneI_act3*ppS2_ppS3_ppS3 + 13.17*geneI_turn)/std::pow(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[70]
    dflux_v50_v_49_dgeneI_inh3 = -ppS2_ppS3_ppS3*(13.17*geneI_act1*ppS2_ppS3_S4 + 13.17*geneI_act2*ppS2_S4_S4 + 13.17*geneI_act3*ppS2_ppS3_ppS3 + 13.17*geneI_turn)/std::pow(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[71]
    dflux_v50_v_49_dgeneI_turn = 13.17/(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[72]
    dflux_v51_v_50_dgeneI_turn = 13.17*geneI;  // dwdp[73]
    dflux_v52_v_51_dgeneJ_inh1 = -ppS2_ppS3_S4*(13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/std::pow(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[74]
    dflux_v52_v_51_dgeneJ_inh2 = -ppS2_S4_S4*(13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/std::pow(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[75]
    dflux_v52_v_51_dgeneJ_inh3 = -ppS2_ppS3_ppS3*(13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/std::pow(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[76]
    dflux_v52_v_51_dgeneJ_turn = 13.17/(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[77]
    dflux_v53_v_52_dgeneJ_turn = 13.17*geneJ;  // dwdp[78]
    dflux_v54_v_53_dgeneK_act1 = 13.17*ppS2_ppS3_S4/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[79]
    dflux_v54_v_53_dgeneK_act2 = 13.17*ppS2_S4_S4/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[80]
    dflux_v54_v_53_dgeneK_inh2 = -ppS2_S4_S4*(13.17*geneK_act1*ppS2_ppS3_S4 + 13.17*geneK_act2*ppS2_S4_S4 + 13.17*geneK_act3*ppS2_ppS3_ppS3 + 13.17*geneK_turn)/std::pow(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[81]
    dflux_v54_v_53_dgeneK_inh3 = -ppS2_ppS3_ppS3*(13.17*geneK_act1*ppS2_ppS3_S4 + 13.17*geneK_act2*ppS2_S4_S4 + 13.17*geneK_act3*ppS2_ppS3_ppS3 + 13.17*geneK_turn)/std::pow(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[82]
    dflux_v54_v_53_dgeneK_turn = 13.17/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[83]
    dflux_v55_v_54_dgeneK_turn = 13.17*geneK;  // dwdp[84]
    dflux_v56_v_55_dgeneL_act2 = 13.17*ppS2_S4_S4/(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[85]
    dflux_v56_v_55_dgeneL_inh1 = -ppS2_ppS3_S4*(13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/std::pow(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[86]
    dflux_v56_v_55_dgeneL_inh2 = -ppS2_S4_S4*(13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/std::pow(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[87]
    dflux_v56_v_55_dgeneL_inh3 = -ppS2_ppS3_ppS3*(13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/std::pow(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdp[88]
    dflux_v56_v_55_dgeneL_turn = 13.17/(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1);  // dwdp[89]
    dflux_v57_v_56_dgeneL_turn = 13.17*geneL;  // dwdp[90]
    dflux_v22_v_21_dk_233 = 13.17*ppS2*std::pow(ppS3, 2);  // dwdp[91]
    dflux_v31_v_30_dk_234 = 13.17*S4*ppS2*ppS3;  // dwdp[92]
    dflux_v27_v_26_dk_244 = 13.17*std::pow(S4, 2)*ppS2;  // dwdp[93]
    dflux_v2_v_1_dpRec_degind = 13.17*TGFb_pRec;  // dwdp[94]
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
