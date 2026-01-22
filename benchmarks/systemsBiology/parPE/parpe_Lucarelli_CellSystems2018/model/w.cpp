#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void w_Lucarelli_CellSystems2018(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    flux_v1_v_0 = 13.17*Rec*Rec_act*TGFb;  // w[0]
    flux_v2_v_1 = 13.17*TGFb_pRec*pRec_degind;  // w[1]
    flux_v3_v_2 = 13.17*S2*std::pow(S4, 2)*k_on_u;  // w[2]
    flux_v4_v_3 = 13.17*S2_S4_S4*kdiss_SS;  // w[3]
    flux_v5_v_4 = 13.17*khomo2*std::pow(ppS2, 3);  // w[4]
    flux_v6_v_5 = 39.509999999999998*S_dephosphos*ppS2_ppS2_ppS2;  // w[5]
    flux_v7_v_6 = 13.17*khomo3*std::pow(ppS3, 3);  // w[6]
    flux_v8_v_7 = 39.509999999999998*S_dephosphos*ppS3_ppS3_ppS3;  // w[7]
    flux_v9_v_8 = 13.17*std::pow(S4, 3)*khomo4;  // w[8]
    flux_v10_v_9 = 13.17*S4_S4_S4*kdiss_SS;  // w[9]
    flux_v11_v_10 = 13.17*S2*S_phos*TGFb_pRec;  // w[10]
    flux_v12_v_11 = 13.17*S_dephosphos*ppS2;  // w[11]
    flux_v13_v_12 = 13.17*S_dephos*pS2;  // w[12]
    flux_v14_v_13 = 13.17*S3*S_phos*TGFb_pRec;  // w[13]
    flux_v15_v_14 = 13.17*S_dephosphos*ppS3;  // w[14]
    flux_v16_v_15 = 13.17*S_dephos*pS3;  // w[15]
    flux_v17_v_16 = 13.17*k_223*std::pow(ppS2, 2)*ppS3;  // w[16]
    flux_v18_v_17 = 26.34*S_dephosphos*ppS2_ppS2_ppS3;  // w[17]
    flux_v19_v_18 = 13.17*S_dephosphos*ppS2_ppS2_ppS3;  // w[18]
    flux_v20_v_19 = 13.17*S4*k_224*std::pow(ppS2, 2);  // w[19]
    flux_v21_v_20 = 26.34*S_dephosphos*ppS2_ppS2_S4;  // w[20]
    flux_v22_v_21 = 13.17*k_233*ppS2*std::pow(ppS3, 2);  // w[21]
    flux_v23_v_22 = 13.17*S_dephosphos*ppS2_ppS3_ppS3;  // w[22]
    flux_v24_v_23 = 26.34*S_dephosphos*ppS2_ppS3_ppS3;  // w[23]
    flux_v25_v_24 = 13.17*S4*k_334*std::pow(ppS3, 2);  // w[24]
    flux_v26_v_25 = 26.34*S_dephosphos*ppS3_ppS3_S4;  // w[25]
    flux_v27_v_26 = 13.17*std::pow(S4, 2)*k_244*ppS2;  // w[26]
    flux_v28_v_27 = 13.17*S_dephosphos*ppS2_S4_S4;  // w[27]
    flux_v29_v_28 = 13.17*std::pow(S4, 2)*k_344*ppS3;  // w[28]
    flux_v30_v_29 = 13.17*S_dephosphos*ppS3_S4_S4;  // w[29]
    flux_v31_v_30 = 13.17*S4*k_234*ppS2*ppS3;  // w[30]
    flux_v32_v_31 = 13.17*S_dephosphos*ppS2_ppS3_S4;  // w[31]
    flux_v33_v_32 = 13.17*S_dephosphos*ppS2_ppS3_S4;  // w[32]
    flux_v34_v_33 = (13.17*geneA_act1*ppS2_ppS3_S4 + 13.17*geneA_act2*ppS2_S4_S4 + 13.17*geneA_act3*ppS2_ppS3_ppS3 + 13.17*geneA_turn)/(geneA_inh1*ppS2_ppS3_S4 + geneA_inh2*ppS2_S4_S4 + geneA_inh3*ppS2_ppS3_ppS3 + 1);  // w[33]
    flux_v35_v_34 = 13.17*geneA*geneA_turn;  // w[34]
    flux_v36_v_35 = (13.17*geneB_act1*ppS2_ppS3_S4 + 13.17*geneB_act2*ppS2_S4_S4 + 13.17*geneB_act3*ppS2_ppS3_ppS3 + 13.17*geneB_turn)/(geneB_inh1*ppS2_ppS3_S4 + geneB_inh2*ppS2_S4_S4 + geneB_inh3*ppS2_ppS3_ppS3 + 1);  // w[35]
    flux_v37_v_36 = 13.17*geneB*geneB_turn;  // w[36]
    flux_v38_v_37 = (13.17*geneC_act1*ppS2_ppS3_S4 + 13.17*geneC_act2*ppS2_S4_S4 + 13.17*geneC_act3*ppS2_ppS3_ppS3 + 13.17*geneC_turn)/(geneC_inh1*ppS2_ppS3_S4 + geneC_inh2*ppS2_S4_S4 + geneC_inh3*ppS2_ppS3_ppS3 + 1);  // w[37]
    flux_v39_v_38 = 13.17*geneC*geneC_turn;  // w[38]
    flux_v40_v_39 = (13.17*geneD_act1*ppS2_ppS3_S4 + 13.17*geneD_act2*ppS2_S4_S4 + 13.17*geneD_act3*ppS2_ppS3_ppS3 + 13.17*geneD_turn)/(geneD_inh1*ppS2_ppS3_S4 + geneD_inh2*ppS2_S4_S4 + geneD_inh3*ppS2_ppS3_ppS3 + 1);  // w[39]
    flux_v41_v_40 = 13.17*geneD*geneD_turn;  // w[40]
    flux_v42_v_41 = (13.17*geneE_act1*ppS2_ppS3_S4 + 13.17*geneE_act2*ppS2_S4_S4 + 13.17*geneE_act3*ppS2_ppS3_ppS3 + 13.17*geneE_turn)/(geneE_inh1*ppS2_ppS3_S4 + geneE_inh2*ppS2_S4_S4 + geneE_inh3*ppS2_ppS3_ppS3 + 1);  // w[41]
    flux_v43_v_42 = 13.17*geneE*geneE_turn;  // w[42]
    flux_v44_v_43 = (13.17*geneF_act1*ppS2_ppS3_S4 + 13.17*geneF_act2*ppS2_S4_S4 + 13.17*geneF_act3*ppS2_ppS3_ppS3 + 13.17*geneF_turn)/(geneF_inh1*ppS2_ppS3_S4 + geneF_inh2*ppS2_S4_S4 + geneF_inh3*ppS2_ppS3_ppS3 + 1);  // w[43]
    flux_v45_v_44 = 13.17*geneF*geneF_turn;  // w[44]
    flux_v46_v_45 = (13.17*geneG_act1*ppS2_ppS3_S4 + 13.17*geneG_act2*ppS2_S4_S4 + 13.17*geneG_act3*ppS2_ppS3_ppS3 + 13.17*geneG_turn)/(geneG_inh1*ppS2_ppS3_S4 + geneG_inh2*ppS2_S4_S4 + geneG_inh3*ppS2_ppS3_ppS3 + 1);  // w[45]
    flux_v47_v_46 = 13.17*geneG*geneG_turn;  // w[46]
    flux_v48_v_47 = (13.17*geneH_act1*ppS2_ppS3_S4 + 13.17*geneH_act2*ppS2_S4_S4 + 13.17*geneH_act3*ppS2_ppS3_ppS3 + 13.17*geneH_turn)/(geneH_inh1*ppS2_ppS3_S4 + geneH_inh2*ppS2_S4_S4 + geneH_inh3*ppS2_ppS3_ppS3 + 1);  // w[47]
    flux_v49_v_48 = 13.17*geneH*geneH_turn;  // w[48]
    flux_v50_v_49 = (13.17*geneI_act1*ppS2_ppS3_S4 + 13.17*geneI_act2*ppS2_S4_S4 + 13.17*geneI_act3*ppS2_ppS3_ppS3 + 13.17*geneI_turn)/(geneI_inh1*ppS2_ppS3_S4 + geneI_inh2*ppS2_S4_S4 + geneI_inh3*ppS2_ppS3_ppS3 + 1);  // w[49]
    flux_v51_v_50 = 13.17*geneI*geneI_turn;  // w[50]
    flux_v52_v_51 = (13.17*geneJ_act1*ppS2_ppS3_S4 + 13.17*geneJ_act2*ppS2_S4_S4 + 13.17*geneJ_act3*ppS2_ppS3_ppS3 + 13.17*geneJ_turn)/(geneJ_inh1*ppS2_ppS3_S4 + geneJ_inh2*ppS2_S4_S4 + geneJ_inh3*ppS2_ppS3_ppS3 + 1);  // w[51]
    flux_v53_v_52 = 13.17*geneJ*geneJ_turn;  // w[52]
    flux_v54_v_53 = (13.17*geneK_act1*ppS2_ppS3_S4 + 13.17*geneK_act2*ppS2_S4_S4 + 13.17*geneK_act3*ppS2_ppS3_ppS3 + 13.17*geneK_turn)/(geneK_inh1*ppS2_ppS3_S4 + geneK_inh2*ppS2_S4_S4 + geneK_inh3*ppS2_ppS3_ppS3 + 1);  // w[53]
    flux_v55_v_54 = 13.17*geneK*geneK_turn;  // w[54]
    flux_v56_v_55 = (13.17*geneL_act1*ppS2_ppS3_S4 + 13.17*geneL_act2*ppS2_S4_S4 + 13.17*geneL_act3*ppS2_ppS3_ppS3 + 13.17*geneL_turn)/(geneL_inh1*ppS2_ppS3_S4 + geneL_inh2*ppS2_S4_S4 + geneL_inh3*ppS2_ppS3_ppS3 + 1);  // w[55]
    flux_v57_v_56 = 13.17*geneL*geneL_turn;  // w[56]
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
