#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Weber_BMC2015 {

// clang-format off

std::array<const char*, 39> parameterNames = {
    "a11", // p[0]
"a12", // p[1]
"a21", // p[2]
"a22", // p[3]
"a31", // p[4]
"a32", // p[5]
"a33", // p[6]
"m11", // p[7]
"m22", // p[8]
"m31", // p[9]
"m33", // p[10]
"p11", // p[11]
"p12", // p[12]
"p13", // p[13]
"p21", // p[14]
"p22", // p[15]
"p31", // p[16]
"p32", // p[17]
"p33", // p[18]
"pu3", // p[19]
"pu4", // p[20]
"pu5", // p[21]
"pu6", // p[22]
"s12", // p[23]
"s21", // p[24]
"s31", // p[25]
"observableParameter1_yCERTpRN24", // p[26]
"observableParameter1_yPI4K3BpRN24", // p[27]
"observableParameter1_yPKDpN0", // p[28]
"observableParameter1_yPKDpN24", // p[29]
"observableParameter1_yPKDpN25", // p[30]
"noiseParameter1_yCERTpRN24", // p[31]
"noiseParameter1_yCERTt", // p[32]
"noiseParameter1_yPI4K3BpRN24", // p[33]
"noiseParameter1_yPI4K3Bt", // p[34]
"noiseParameter1_yPKDpN0", // p[35]
"noiseParameter1_yPKDpN24", // p[36]
"noiseParameter1_yPKDpN25", // p[37]
"noiseParameter1_yPKDt", // p[38]
};

std::array<const char*, 7> fixedParameterNames = {
    "pu2", // k[0]
"Ect_Expr_CERT_flag", // k[1]
"Ect_Expr_PI4K3beta_flag", // k[2]
"PdBu_dose", // k[3]
"PdBu_time", // k[4]
"kb_NB142_70_dose", // k[5]
"kb_NB142_70_time", // k[6]
};

std::array<const char*, 7> stateNames = {
    "PKD", // x_rdata[0]
"PKD_p:DAG", // x_rdata[1]
"PI4KIII\beta", // x_rdata[2]
"PI4KIII\beta_p", // x_rdata[3]
"CERT_{a,ER}", // x_rdata[4]
"CERT_{p,ER}", // x_rdata[5]
"CERT_{a,TGN}", // x_rdata[6]
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

std::array<const char*, 26> expressionNames = {
    "u2", // w[0]
"u3", // w[1]
"u4", // w[2]
"u5", // w[3]
"u6", // w[4]
"flux_v1_v_0", // w[5]
"flux_v2_v_1", // w[6]
"flux_v3_v_2", // w[7]
"flux_v4_v_3", // w[8]
"flux_v5_v_4", // w[9]
"flux_v6_v_5", // w[10]
"flux_v7_v_6", // w[11]
"flux_v8_v_7", // w[12]
"flux_v9_v_8", // w[13]
"flux_v10_v_9", // w[14]
"flux_v11_v_10", // w[15]
"flux_v12_v_11", // w[16]
"flux_v13_v_12", // w[17]
"flux_v14_v_13", // w[18]
"flux_v15_v_14", // w[19]
"flux_v16_v_15", // w[20]
"flux_v17_v_16", // w[21]
"flux_v18_v_17", // w[22]
"flux_v19_v_18", // w[23]
"flux_v20_v_19", // w[24]
"flux_v21_v_20", // w[25]
};

std::array<const char*, 39> parameterIds = {
    "a11", // p[0]
"a12", // p[1]
"a21", // p[2]
"a22", // p[3]
"a31", // p[4]
"a32", // p[5]
"a33", // p[6]
"m11", // p[7]
"m22", // p[8]
"m31", // p[9]
"m33", // p[10]
"p11", // p[11]
"p12", // p[12]
"p13", // p[13]
"p21", // p[14]
"p22", // p[15]
"p31", // p[16]
"p32", // p[17]
"p33", // p[18]
"pu3", // p[19]
"pu4", // p[20]
"pu5", // p[21]
"pu6", // p[22]
"s12", // p[23]
"s21", // p[24]
"s31", // p[25]
"observableParameter1_yCERTpRN24", // p[26]
"observableParameter1_yPI4K3BpRN24", // p[27]
"observableParameter1_yPKDpN0", // p[28]
"observableParameter1_yPKDpN24", // p[29]
"observableParameter1_yPKDpN25", // p[30]
"noiseParameter1_yCERTpRN24", // p[31]
"noiseParameter1_yCERTt", // p[32]
"noiseParameter1_yPI4K3BpRN24", // p[33]
"noiseParameter1_yPI4K3Bt", // p[34]
"noiseParameter1_yPKDpN0", // p[35]
"noiseParameter1_yPKDpN24", // p[36]
"noiseParameter1_yPKDpN25", // p[37]
"noiseParameter1_yPKDt", // p[38]
};

std::array<const char*, 7> fixedParameterIds = {
    "pu2", // k[0]
"Ect_Expr_CERT_flag", // k[1]
"Ect_Expr_PI4K3beta_flag", // k[2]
"PdBu_dose", // k[3]
"PdBu_time", // k[4]
"kb_NB142_70_dose", // k[5]
"kb_NB142_70_time", // k[6]
};

std::array<const char*, 7> stateIds = {
    "PKD", // x_rdata[0]
"PKDDAGa", // x_rdata[1]
"PI4K3B", // x_rdata[2]
"PI4K3Ba", // x_rdata[3]
"CERTERa", // x_rdata[4]
"CERT", // x_rdata[5]
"CERTTGNa", // x_rdata[6]
};

std::array<const char*, 8> observableIds = {
    "yCERTpRN24", // y[0]
"yCERTt", // y[1]
"yPI4K3BpRN24", // y[2]
"yPI4K3Bt", // y[3]
"yPKDpN0", // y[4]
"yPKDpN24", // y[5]
"yPKDpN25", // y[6]
"yPKDt", // y[7]
};

std::array<const char*, 26> expressionIds = {
    "u2", // w[0]
"u3", // w[1]
"u4", // w[2]
"u5", // w[3]
"u6", // w[4]
"flux_v1_v_0", // w[5]
"flux_v2_v_1", // w[6]
"flux_v3_v_2", // w[7]
"flux_v4_v_3", // w[8]
"flux_v5_v_4", // w[9]
"flux_v6_v_5", // w[10]
"flux_v7_v_6", // w[11]
"flux_v8_v_7", // w[12]
"flux_v9_v_8", // w[13]
"flux_v10_v_9", // w[14]
"flux_v11_v_10", // w[15]
"flux_v12_v_11", // w[16]
"flux_v13_v_12", // w[17]
"flux_v14_v_13", // w[18]
"flux_v15_v_14", // w[19]
"flux_v16_v_15", // w[20]
"flux_v17_v_16", // w[21]
"flux_v18_v_17", // w[22]
"flux_v19_v_18", // w[23]
"flux_v20_v_19", // w[24]
"flux_v21_v_20", // w[25]
};

std::array<int, 7> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6
};

std::array<bool, 6> rootInitialValues = {
    true, true, true, true, true, true
};

// clang-format on

} // namespace model_Weber_BMC2015

} // namespace amici
