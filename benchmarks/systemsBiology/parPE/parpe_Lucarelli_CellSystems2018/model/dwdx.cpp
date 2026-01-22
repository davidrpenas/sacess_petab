#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<sunindextype, 34> dwdx_colptrs_Lucarelli_CellSystems2018_ = {
    0, 1, 2, 5, 7, 8, 15, 16, 17, 18, 19, 20, 21, 28, 35, 36, 38, 52, 53, 67, 68, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93
};

void dwdx_colptrs_Lucarelli_CellSystems2018(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Lucarelli_CellSystems2018_));
}
} // namespace model_Lucarelli_CellSystems2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<sunindextype, 93> dwdx_rowvals_Lucarelli_CellSystems2018_ = {
    0, 0, 1, 10, 13, 2, 10, 13, 2, 8, 19, 24, 26, 28, 30, 3, 5, 7, 9, 12, 15, 4, 11, 16, 19, 21, 26, 30, 6, 14, 16, 21, 24, 28, 30, 20, 17, 18, 22, 23, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 25, 31, 32, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 29, 27, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56
};

void dwdx_rowvals_Lucarelli_CellSystems2018(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Lucarelli_CellSystems2018_));
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
#include "dwdx.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void dwdx_Lucarelli_CellSystems2018(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v2_v_1_dTGFb_pRec = 13.17*pRec_degind;  // dwdx[2]
        dflux_v4_v_3_dS2_S4_S4 = 13.17*kdiss_SS;  // dwdx[15]
        dflux_v6_v_5_dppS2_ppS2_ppS2 = 39.509999999999998*S_dephosphos;  // dwdx[16]
        dflux_v8_v_7_dppS3_ppS3_ppS3 = 39.509999999999998*S_dephosphos;  // dwdx[17]
        dflux_v10_v_9_dS4_S4_S4 = 13.17*kdiss_SS;  // dwdx[18]
        dflux_v13_v_12_dpS2 = 13.17*S_dephos;  // dwdx[19]
        dflux_v16_v_15_dpS3 = 13.17*S_dephos;  // dwdx[20]
        dflux_v12_v_11_dppS2 = 13.17*S_dephosphos;  // dwdx[22]
        dflux_v15_v_14_dppS3 = 13.17*S_dephosphos;  // dwdx[29]
        dflux_v21_v_20_dppS2_ppS2_S4 = 26.34*S_dephosphos;  // dwdx[35]
        dflux_v18_v_17_dppS2_ppS2_ppS3 = 26.34*S_dephosphos;  // dwdx[36]
        dflux_v19_v_18_dppS2_ppS2_ppS3 = 13.17*S_dephosphos;  // dwdx[37]
        dflux_v23_v_22_dppS2_ppS3_ppS3 = 13.17*S_dephosphos;  // dwdx[38]
        dflux_v24_v_23_dppS2_ppS3_ppS3 = 26.34*S_dephosphos;  // dwdx[39]
        dflux_v26_v_25_dppS3_ppS3_S4 = 26.34*S_dephosphos;  // dwdx[52]
        dflux_v32_v_31_dppS2_ppS3_S4 = 13.17*S_dephosphos;  // dwdx[53]
        dflux_v33_v_32_dppS2_ppS3_S4 = 13.17*S_dephosphos;  // dwdx[54]
        dflux_v30_v_29_dppS3_S4_S4 = 13.17*S_dephosphos;  // dwdx[67]
        dflux_v28_v_27_dppS2_S4_S4 = 13.17*S_dephosphos;  // dwdx[68]
        dflux_v35_v_34_dgeneA = 13.17*geneA_turn;  // dwdx[81]
        dflux_v37_v_36_dgeneB = 13.17*geneB_turn;  // dwdx[82]
        dflux_v39_v_38_dgeneC = 13.17*geneC_turn;  // dwdx[83]
        dflux_v41_v_40_dgeneD = 13.17*geneD_turn;  // dwdx[84]
        dflux_v43_v_42_dgeneE = 13.17*geneE_turn;  // dwdx[85]
        dflux_v45_v_44_dgeneF = 13.17*geneF_turn;  // dwdx[86]
        dflux_v47_v_46_dgeneG = 13.17*geneG_turn;  // dwdx[87]
        dflux_v49_v_48_dgeneH = 13.17*geneH_turn;  // dwdx[88]
        dflux_v51_v_50_dgeneI = 13.17*geneI_turn;  // dwdx[89]
        dflux_v53_v_52_dgeneJ = 13.17*geneJ_turn;  // dwdx[90]
        dflux_v55_v_54_dgeneK = 13.17*geneK_turn;  // dwdx[91]
        dflux_v57_v_56_dgeneL = 13.17*geneL_turn;  // dwdx[92]
    }

    // dynamic expressions
    dflux_v1_v_0_dTGFb = 13.17*Rec*Rec_act;  // dwdx[0]
    dflux_v1_v_0_dRec = 13.17*Rec_act*TGFb;  // dwdx[1]
    dflux_v11_v_10_dTGFb_pRec = 13.17*S2*S_phos;  // dwdx[3]
    dflux_v14_v_13_dTGFb_pRec = 13.17*S3*S_phos;  // dwdx[4]
    dflux_v3_v_2_dS2 = 13.17*std::pow(S4, 2)*k_on_u;  // dwdx[5]
    dflux_v11_v_10_dS2 = 13.17*S_phos*TGFb_pRec;  // dwdx[6]
    dflux_v14_v_13_dS3 = 13.17*S_phos*TGFb_pRec;  // dwdx[7]
    dflux_v3_v_2_dS4 = 26.34*S2*S4*k_on_u;  // dwdx[8]
    dflux_v9_v_8_dS4 = 39.509999999999998*std::pow(S4, 2)*khomo4;  // dwdx[9]
    dflux_v20_v_19_dS4 = 13.17*k_224*std::pow(ppS2, 2);  // dwdx[10]
    dflux_v25_v_24_dS4 = 13.17*k_334*std::pow(ppS3, 2);  // dwdx[11]
    dflux_v27_v_26_dS4 = 26.34*S4*k_244*ppS2;  // dwdx[12]
    dflux_v29_v_28_dS4 = 26.34*S4*k_344*ppS3;  // dwdx[13]
    dflux_v31_v_30_dS4 = 13.17*k_234*ppS2*ppS3;  // dwdx[14]
    dflux_v5_v_4_dppS2 = 39.509999999999998*khomo2*std::pow(ppS2, 2);  // dwdx[21]
    dflux_v17_v_16_dppS2 = 26.34*k_223*ppS2*ppS3;  // dwdx[23]
    dflux_v20_v_19_dppS2 = 26.34*S4*k_224*ppS2;  // dwdx[24]
    dflux_v22_v_21_dppS2 = 13.17*k_233*std::pow(ppS3, 2);  // dwdx[25]
    dflux_v27_v_26_dppS2 = 13.17*std::pow(S4, 2)*k_244;  // dwdx[26]
    dflux_v31_v_30_dppS2 = 13.17*S4*k_234*ppS3;  // dwdx[27]
    dflux_v7_v_6_dppS3 = 39.509999999999998*khomo3*std::pow(ppS3, 2);  // dwdx[28]
    dflux_v17_v_16_dppS3 = 13.17*k_223*std::pow(ppS2, 2);  // dwdx[30]
    dflux_v22_v_21_dppS3 = 26.34*k_233*ppS2*ppS3;  // dwdx[31]
    dflux_v25_v_24_dppS3 = 26.34*S4*k_334*ppS3;  // dwdx[32]
    dflux_v29_v_28_dppS3 = 13.17*std::pow(S4, 2)*k_344;  // dwdx[33]
    dflux_v31_v_30_dppS3 = 13.17*S4*k_234*ppS2;  // dwdx[34]
    dflux_v34_v_33_dppS2_ppS3_ppS3 = 13.17*geneA_act3/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1) - geneA_inh3*(13.17*geneA_act1*ppS2_ppS3_S4 + 13.17*geneA_act2*ppS2_S4_S4 + 13.17*geneA_act3*ppS2_ppS3_ppS3 + 13.17*geneA_turn)/std::pow(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[40]
    dflux_v36_v_35_dppS2_ppS3_ppS3 = 13.17*geneB_act3/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1) - geneB_inh3*(13.17*geneB_act1*ppS2_ppS3_S4 + 13.17*geneB_act2*ppS2_S4_S4 + 13.17*geneB_act3*ppS2_ppS3_ppS3 + 13.17*geneB_turn)/std::pow(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[41]
    dflux_v38_v_37_dppS2_ppS3_ppS3 = 13.17*geneC_act3/(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1) - geneC_inh3*(13.17*geneC_act1*ppS2_ppS3_S4 + 13.17*geneC_act2*ppS2_S4_S4 + 13.17*geneC_act3*ppS2_ppS3_ppS3 + 13.17*geneC_turn)/std::pow(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[42]
    dflux_v40_v_39_dppS2_ppS3_ppS3 = 13.17*geneD_act3/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1) - geneD_inh3*(13.17*geneD_act1*ppS2_ppS3_S4 + 13.17*geneD_act2*ppS2_S4_S4 + 13.17*geneD_act3*ppS2_ppS3_ppS3 + 13.17*geneD_turn)/std::pow(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[43]
    dflux_v42_v_41_dppS2_ppS3_ppS3 = 13.17*geneE_act3/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1) - geneE_inh3*(13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/std::pow(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[44]
    dflux_v44_v_43_dppS2_ppS3_ppS3 = 13.17*geneF_act3/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1) - geneF_inh3*(13.17*geneF_act1*ppS2_ppS3_S4 + 13.17*geneF_act2*ppS2_S4_S4 + 13.17*geneF_act3*ppS2_ppS3_ppS3 + 13.17*geneF_turn)/std::pow(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[45]
    dflux_v46_v_45_dppS2_ppS3_ppS3 = 13.17*geneG_act3/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1) - geneG_inh3*(13.17*geneG_act1*ppS2_ppS3_S4 + 13.17*geneG_act2*ppS2_S4_S4 + 13.17*geneG_act3*ppS2_ppS3_ppS3 + 13.17*geneG_turn)/std::pow(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[46]
    dflux_v48_v_47_dppS2_ppS3_ppS3 = 13.17*geneH_act3/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1) - geneH_inh3*(13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/std::pow(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[47]
    dflux_v50_v_49_dppS2_ppS3_ppS3 = 13.17*geneI_act3/(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1) - geneI_inh3*(13.17*geneI_act1*ppS2_ppS3_S4 + 13.17*geneI_act2*ppS2_S4_S4 + 13.17*geneI_act3*ppS2_ppS3_ppS3 + 13.17*geneI_turn)/std::pow(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[48]
    dflux_v52_v_51_dppS2_ppS3_ppS3 = 13.17*geneJ_act3/(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1) - geneJ_inh3*(13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/std::pow(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[49]
    dflux_v54_v_53_dppS2_ppS3_ppS3 = 13.17*geneK_act3/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1) - geneK_inh3*(13.17*geneK_act1*ppS2_ppS3_S4 + 13.17*geneK_act2*ppS2_S4_S4 + 13.17*geneK_act3*ppS2_ppS3_ppS3 + 13.17*geneK_turn)/std::pow(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[50]
    dflux_v56_v_55_dppS2_ppS3_ppS3 = 13.17*geneL_act3/(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1) - geneL_inh3*(13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/std::pow(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[51]
    dflux_v34_v_33_dppS2_ppS3_S4 = 13.17*geneA_act1/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1) - geneA_inh1*(13.17*geneA_act1*ppS2_ppS3_S4 + 13.17*geneA_act2*ppS2_S4_S4 + 13.17*geneA_act3*ppS2_ppS3_ppS3 + 13.17*geneA_turn)/std::pow(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[55]
    dflux_v36_v_35_dppS2_ppS3_S4 = 13.17*geneB_act1/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1) - geneB_inh1*(13.17*geneB_act1*ppS2_ppS3_S4 + 13.17*geneB_act2*ppS2_S4_S4 + 13.17*geneB_act3*ppS2_ppS3_ppS3 + 13.17*geneB_turn)/std::pow(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[56]
    dflux_v38_v_37_dppS2_ppS3_S4 = 13.17*geneC_act1/(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1) - geneC_inh1*(13.17*geneC_act1*ppS2_ppS3_S4 + 13.17*geneC_act2*ppS2_S4_S4 + 13.17*geneC_act3*ppS2_ppS3_ppS3 + 13.17*geneC_turn)/std::pow(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[57]
    dflux_v40_v_39_dppS2_ppS3_S4 = 13.17*geneD_act1/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1) - geneD_inh1*(13.17*geneD_act1*ppS2_ppS3_S4 + 13.17*geneD_act2*ppS2_S4_S4 + 13.17*geneD_act3*ppS2_ppS3_ppS3 + 13.17*geneD_turn)/std::pow(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[58]
    dflux_v42_v_41_dppS2_ppS3_S4 = 13.17*geneE_act1/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1) - geneE_inh1*(13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/std::pow(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[59]
    dflux_v44_v_43_dppS2_ppS3_S4 = 13.17*geneF_act1/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1) - geneF_inh1*(13.17*geneF_act1*ppS2_ppS3_S4 + 13.17*geneF_act2*ppS2_S4_S4 + 13.17*geneF_act3*ppS2_ppS3_ppS3 + 13.17*geneF_turn)/std::pow(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[60]
    dflux_v46_v_45_dppS2_ppS3_S4 = 13.17*geneG_act1/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1) - geneG_inh1*(13.17*geneG_act1*ppS2_ppS3_S4 + 13.17*geneG_act2*ppS2_S4_S4 + 13.17*geneG_act3*ppS2_ppS3_ppS3 + 13.17*geneG_turn)/std::pow(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[61]
    dflux_v48_v_47_dppS2_ppS3_S4 = 13.17*geneH_act1/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1) - geneH_inh1*(13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/std::pow(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[62]
    dflux_v50_v_49_dppS2_ppS3_S4 = 13.17*geneI_act1/(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1) - geneI_inh1*(13.17*geneI_act1*ppS2_ppS3_S4 + 13.17*geneI_act2*ppS2_S4_S4 + 13.17*geneI_act3*ppS2_ppS3_ppS3 + 13.17*geneI_turn)/std::pow(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[63]
    dflux_v52_v_51_dppS2_ppS3_S4 = 13.17*geneJ_act1/(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1) - geneJ_inh1*(13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/std::pow(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[64]
    dflux_v54_v_53_dppS2_ppS3_S4 = 13.17*geneK_act1/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1) - geneK_inh1*(13.17*geneK_act1*ppS2_ppS3_S4 + 13.17*geneK_act2*ppS2_S4_S4 + 13.17*geneK_act3*ppS2_ppS3_ppS3 + 13.17*geneK_turn)/std::pow(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[65]
    dflux_v56_v_55_dppS2_ppS3_S4 = 13.17*geneL_act1/(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1) - geneL_inh1*(13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/std::pow(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[66]
    dflux_v34_v_33_dppS2_S4_S4 = 13.17*geneA_act2/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1) - geneA_inh2*(13.17*geneA_act1*ppS2_ppS3_S4 + 13.17*geneA_act2*ppS2_S4_S4 + 13.17*geneA_act3*ppS2_ppS3_ppS3 + 13.17*geneA_turn)/std::pow(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[69]
    dflux_v36_v_35_dppS2_S4_S4 = 13.17*geneB_act2/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1) - geneB_inh2*(13.17*geneB_act1*ppS2_ppS3_S4 + 13.17*geneB_act2*ppS2_S4_S4 + 13.17*geneB_act3*ppS2_ppS3_ppS3 + 13.17*geneB_turn)/std::pow(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[70]
    dflux_v38_v_37_dppS2_S4_S4 = 13.17*geneC_act2/(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1) - geneC_inh2*(13.17*geneC_act1*ppS2_ppS3_S4 + 13.17*geneC_act2*ppS2_S4_S4 + 13.17*geneC_act3*ppS2_ppS3_ppS3 + 13.17*geneC_turn)/std::pow(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[71]
    dflux_v40_v_39_dppS2_S4_S4 = 13.17*geneD_act2/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1) - geneD_inh2*(13.17*geneD_act1*ppS2_ppS3_S4 + 13.17*geneD_act2*ppS2_S4_S4 + 13.17*geneD_act3*ppS2_ppS3_ppS3 + 13.17*geneD_turn)/std::pow(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[72]
    dflux_v42_v_41_dppS2_S4_S4 = 13.17*geneE_act2/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1) - geneE_inh2*(13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/std::pow(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[73]
    dflux_v44_v_43_dppS2_S4_S4 = 13.17*geneF_act2/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1) - geneF_inh2*(13.17*geneF_act1*ppS2_ppS3_S4 + 13.17*geneF_act2*ppS2_S4_S4 + 13.17*geneF_act3*ppS2_ppS3_ppS3 + 13.17*geneF_turn)/std::pow(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[74]
    dflux_v46_v_45_dppS2_S4_S4 = 13.17*geneG_act2/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1) - geneG_inh2*(13.17*geneG_act1*ppS2_ppS3_S4 + 13.17*geneG_act2*ppS2_S4_S4 + 13.17*geneG_act3*ppS2_ppS3_ppS3 + 13.17*geneG_turn)/std::pow(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[75]
    dflux_v48_v_47_dppS2_S4_S4 = 13.17*geneH_act2/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1) - geneH_inh2*(13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/std::pow(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[76]
    dflux_v50_v_49_dppS2_S4_S4 = 13.17*geneI_act2/(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1) - geneI_inh2*(13.17*geneI_act1*ppS2_ppS3_S4 + 13.17*geneI_act2*ppS2_S4_S4 + 13.17*geneI_act3*ppS2_ppS3_ppS3 + 13.17*geneI_turn)/std::pow(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[77]
    dflux_v52_v_51_dppS2_S4_S4 = 13.17*geneJ_act2/(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1) - geneJ_inh2*(13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/std::pow(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[78]
    dflux_v54_v_53_dppS2_S4_S4 = 13.17*geneK_act2/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1) - geneK_inh2*(13.17*geneK_act1*ppS2_ppS3_S4 + 13.17*geneK_act2*ppS2_S4_S4 + 13.17*geneK_act3*ppS2_ppS3_ppS3 + 13.17*geneK_turn)/std::pow(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[79]
    dflux_v56_v_55_dppS2_S4_S4 = 13.17*geneL_act2/(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1) - geneL_inh2*(13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/std::pow(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1, 2);  // dwdx[80]
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
