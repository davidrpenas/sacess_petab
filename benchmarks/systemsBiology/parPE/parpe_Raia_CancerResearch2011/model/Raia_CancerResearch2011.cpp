#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Raia_CancerResearch2011 {

// clang-format off

std::array<const char*, 39> parameterNames = {
    "CD274mRNA_production", // p[0]
"DecoyR_binding", // p[1]
"JAK2_p_inhibition", // p[2]
"JAK2_phosphorylation", // p[3]
"Kon_IL13Rec", // p[4]
"Rec_intern", // p[5]
"Rec_phosphorylation", // p[6]
"Rec_recycle", // p[7]
"SOCS3_accumulation", // p[8]
"SOCS3_degradation", // p[9]
"SOCS3_translation", // p[10]
"SOCS3mRNA_production", // p[11]
"STAT5_phosphorylation", // p[12]
"init_Rec_i", // p[13]
"pJAK2_dephosphorylation", // p[14]
"pRec_degradation", // p[15]
"pRec_intern", // p[16]
"pSTAT5_dephosphorylation", // p[17]
"observableParameter1_observable_CD274mRNA", // p[18]
"observableParameter1_observable_IL13_cell", // p[19]
"observableParameter1_observable_SOCS3mRNA", // p[20]
"observableParameter1_observable_pIL4Ra", // p[21]
"observableParameter1_observable_pJAK2", // p[22]
"noiseParameter1_observable_CD274mRNA", // p[23]
"noiseParameter2_observable_CD274mRNA", // p[24]
"noiseParameter1_observable_IL13_cell", // p[25]
"noiseParameter2_observable_IL13_cell", // p[26]
"noiseParameter1_observable_RecSurf", // p[27]
"noiseParameter2_observable_RecSurf", // p[28]
"noiseParameter1_observable_SOCS3", // p[29]
"noiseParameter2_observable_SOCS3", // p[30]
"noiseParameter1_observable_SOCS3mRNA", // p[31]
"noiseParameter2_observable_SOCS3mRNA", // p[32]
"noiseParameter1_observable_pIL4Ra", // p[33]
"noiseParameter2_observable_pIL4Ra", // p[34]
"noiseParameter1_observable_pJAK2", // p[35]
"noiseParameter2_observable_pJAK2", // p[36]
"noiseParameter1_observable_pSTAT5", // p[37]
"noiseParameter2_observable_pSTAT5", // p[38]
};

std::array<const char*, 1> fixedParameterNames = {
    "il13_level", // k[0]
};

std::array<const char*, 14> stateNames = {
    "Rec", // x_rdata[0]
"Rec_i", // x_rdata[1]
"IL13_Rec", // x_rdata[2]
"p_IL13_Rec", // x_rdata[3]
"p_IL13_Rec_i", // x_rdata[4]
"JAK2", // x_rdata[5]
"pJAK2", // x_rdata[6]
"STAT5", // x_rdata[7]
"pSTAT5", // x_rdata[8]
"SOCS3mRNA", // x_rdata[9]
"DecoyR", // x_rdata[10]
"IL13_DecoyR", // x_rdata[11]
"SOCS3", // x_rdata[12]
"CD274mRNA", // x_rdata[13]
};

std::array<const char*, 8> observableNames = {
    "", // y[0]
"", // y[1]
"", // y[2]
"", // y[3]
"", // y[4]
"", // y[5]
"", // y[6]
"", // y[7]
};

std::array<const ObservableScaling, 8> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
ObservableScaling::lin, // y[3]
ObservableScaling::lin, // y[4]
ObservableScaling::lin, // y[5]
ObservableScaling::lin, // y[6]
ObservableScaling::lin, // y[7]
};

std::array<const char*, 16> expressionNames = {
    "flux_v1_v_0", // w[0]
"flux_v2_v_1", // w[1]
"flux_v3_v_2", // w[2]
"flux_v4_v_3", // w[3]
"flux_v5_v_4", // w[4]
"flux_v6_v_5", // w[5]
"flux_v7_v_6", // w[6]
"flux_v8_v_7", // w[7]
"flux_v9_v_8", // w[8]
"flux_v10_v_9", // w[9]
"flux_v11_v_10", // w[10]
"flux_v12_v_11", // w[11]
"flux_v13_v_12", // w[12]
"flux_v14_v_13", // w[13]
"flux_v15_v_14", // w[14]
"flux_v16_v_15", // w[15]
};

std::array<const char*, 39> parameterIds = {
    "CD274mRNA_production", // p[0]
"DecoyR_binding", // p[1]
"JAK2_p_inhibition", // p[2]
"JAK2_phosphorylation", // p[3]
"Kon_IL13Rec", // p[4]
"Rec_intern", // p[5]
"Rec_phosphorylation", // p[6]
"Rec_recycle", // p[7]
"SOCS3_accumulation", // p[8]
"SOCS3_degradation", // p[9]
"SOCS3_translation", // p[10]
"SOCS3mRNA_production", // p[11]
"STAT5_phosphorylation", // p[12]
"init_Rec_i", // p[13]
"pJAK2_dephosphorylation", // p[14]
"pRec_degradation", // p[15]
"pRec_intern", // p[16]
"pSTAT5_dephosphorylation", // p[17]
"observableParameter1_observable_CD274mRNA", // p[18]
"observableParameter1_observable_IL13_cell", // p[19]
"observableParameter1_observable_SOCS3mRNA", // p[20]
"observableParameter1_observable_pIL4Ra", // p[21]
"observableParameter1_observable_pJAK2", // p[22]
"noiseParameter1_observable_CD274mRNA", // p[23]
"noiseParameter2_observable_CD274mRNA", // p[24]
"noiseParameter1_observable_IL13_cell", // p[25]
"noiseParameter2_observable_IL13_cell", // p[26]
"noiseParameter1_observable_RecSurf", // p[27]
"noiseParameter2_observable_RecSurf", // p[28]
"noiseParameter1_observable_SOCS3", // p[29]
"noiseParameter2_observable_SOCS3", // p[30]
"noiseParameter1_observable_SOCS3mRNA", // p[31]
"noiseParameter2_observable_SOCS3mRNA", // p[32]
"noiseParameter1_observable_pIL4Ra", // p[33]
"noiseParameter2_observable_pIL4Ra", // p[34]
"noiseParameter1_observable_pJAK2", // p[35]
"noiseParameter2_observable_pJAK2", // p[36]
"noiseParameter1_observable_pSTAT5", // p[37]
"noiseParameter2_observable_pSTAT5", // p[38]
};

std::array<const char*, 1> fixedParameterIds = {
    "il13_level", // k[0]
};

std::array<const char*, 14> stateIds = {
    "Rec", // x_rdata[0]
"Rec_i", // x_rdata[1]
"IL13_Rec", // x_rdata[2]
"p_IL13_Rec", // x_rdata[3]
"p_IL13_Rec_i", // x_rdata[4]
"JAK2", // x_rdata[5]
"pJAK2", // x_rdata[6]
"STAT5", // x_rdata[7]
"pSTAT5", // x_rdata[8]
"SOCS3mRNA", // x_rdata[9]
"DecoyR", // x_rdata[10]
"IL13_DecoyR", // x_rdata[11]
"SOCS3", // x_rdata[12]
"CD274mRNA", // x_rdata[13]
};

std::array<const char*, 8> observableIds = {
    "observable_CD274mRNA", // y[0]
"observable_IL13_cell", // y[1]
"observable_RecSurf", // y[2]
"observable_SOCS3", // y[3]
"observable_SOCS3mRNA", // y[4]
"observable_pIL4Ra", // y[5]
"observable_pJAK2", // y[6]
"observable_pSTAT5", // y[7]
};

std::array<const char*, 16> expressionIds = {
    "flux_v1_v_0", // w[0]
"flux_v2_v_1", // w[1]
"flux_v3_v_2", // w[2]
"flux_v4_v_3", // w[3]
"flux_v5_v_4", // w[4]
"flux_v6_v_5", // w[5]
"flux_v7_v_6", // w[6]
"flux_v8_v_7", // w[7]
"flux_v9_v_8", // w[8]
"flux_v10_v_9", // w[9]
"flux_v11_v_10", // w[10]
"flux_v12_v_11", // w[11]
"flux_v13_v_12", // w[12]
"flux_v14_v_13", // w[13]
"flux_v15_v_14", // w[14]
"flux_v16_v_15", // w[15]
};

std::array<int, 14> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Raia_CancerResearch2011

} // namespace amici
