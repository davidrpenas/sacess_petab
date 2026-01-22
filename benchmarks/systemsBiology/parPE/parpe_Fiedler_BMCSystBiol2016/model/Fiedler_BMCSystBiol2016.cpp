#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Fiedler_BMCSystBiol2016 {

// clang-format off

std::array<const char*, 48> parameterNames = {
    "K_1", // p[0]
"K_2", // p[1]
"K_3", // p[2]
"k10", // p[3]
"k11", // p[4]
"k2", // p[5]
"k3", // p[6]
"k4", // p[7]
"k5", // p[8]
"k6", // p[9]
"tau1", // p[10]
"tau2", // p[11]
"observableParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // p[12]
"observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // p[13]
"observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // p[14]
"observableParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // p[15]
"observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // p[16]
"observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // p[17]
"observableParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // p[18]
"observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // p[19]
"observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // p[20]
"observableParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // p[21]
"observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // p[22]
"observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // p[23]
"noiseParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // p[24]
"noiseParameter2_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // p[25]
"noiseParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // p[26]
"noiseParameter2_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // p[27]
"noiseParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // p[28]
"noiseParameter2_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // p[29]
"noiseParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // p[30]
"noiseParameter2_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // p[31]
"noiseParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // p[32]
"noiseParameter2_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // p[33]
"noiseParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // p[34]
"noiseParameter2_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // p[35]
"noiseParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // p[36]
"noiseParameter2_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // p[37]
"noiseParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // p[38]
"noiseParameter2_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // p[39]
"noiseParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // p[40]
"noiseParameter2_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // p[41]
"noiseParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // p[42]
"noiseParameter2_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // p[43]
"noiseParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // p[44]
"noiseParameter2_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // p[45]
"noiseParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // p[46]
"noiseParameter2_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // p[47]
};

std::array<const char*, 5> fixedParameterNames = {
    "ERK_total", // k[0]
"MEK_total", // k[1]
"RAF_total", // k[2]
"Sorafenib", // k[3]
"UO126", // k[4]
};

std::array<const char*, 6> stateNames = {
    "RAF", // x_rdata[0]
"pRAF", // x_rdata[1]
"MEK", // x_rdata[2]
"pMEK", // x_rdata[3]
"ERK", // x_rdata[4]
"pERK", // x_rdata[5]
};

std::array<const char*, 12> observableNames = {
    "", // y[0]
"", // y[1]
"", // y[2]
"", // y[3]
"", // y[4]
"", // y[5]
"", // y[6]
"", // y[7]
"", // y[8]
"", // y[9]
"", // y[10]
"", // y[11]
};

std::array<const ObservableScaling, 12> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
ObservableScaling::lin, // y[3]
ObservableScaling::lin, // y[4]
ObservableScaling::lin, // y[5]
ObservableScaling::lin, // y[6]
ObservableScaling::lin, // y[7]
ObservableScaling::lin, // y[8]
ObservableScaling::lin, // y[9]
ObservableScaling::lin, // y[10]
ObservableScaling::lin, // y[11]
};

std::array<const char*, 7> expressionNames = {
    "k1max", // w[0]
"flux_v1_v_0", // w[1]
"flux_v2_v_1", // w[2]
"flux_v3_v_2", // w[3]
"flux_v4_v_3", // w[4]
"flux_v5_v_4", // w[5]
"flux_v6_v_5", // w[6]
};

std::array<const char*, 48> parameterIds = {
    "K_1", // p[0]
"K_2", // p[1]
"K_3", // p[2]
"k10", // p[3]
"k11", // p[4]
"k2", // p[5]
"k3", // p[6]
"k4", // p[7]
"k5", // p[8]
"k6", // p[9]
"tau1", // p[10]
"tau2", // p[11]
"observableParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // p[12]
"observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // p[13]
"observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // p[14]
"observableParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // p[15]
"observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // p[16]
"observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // p[17]
"observableParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // p[18]
"observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // p[19]
"observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // p[20]
"observableParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // p[21]
"observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // p[22]
"observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // p[23]
"noiseParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // p[24]
"noiseParameter2_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // p[25]
"noiseParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // p[26]
"noiseParameter2_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // p[27]
"noiseParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // p[28]
"noiseParameter2_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // p[29]
"noiseParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // p[30]
"noiseParameter2_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // p[31]
"noiseParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // p[32]
"noiseParameter2_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // p[33]
"noiseParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // p[34]
"noiseParameter2_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // p[35]
"noiseParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // p[36]
"noiseParameter2_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // p[37]
"noiseParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // p[38]
"noiseParameter2_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // p[39]
"noiseParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // p[40]
"noiseParameter2_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // p[41]
"noiseParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // p[42]
"noiseParameter2_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // p[43]
"noiseParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // p[44]
"noiseParameter2_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // p[45]
"noiseParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // p[46]
"noiseParameter2_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // p[47]
};

std::array<const char*, 5> fixedParameterIds = {
    "ERK_total", // k[0]
"MEK_total", // k[1]
"RAF_total", // k[2]
"Sorafenib", // k[3]
"UO126", // k[4]
};

std::array<const char*, 6> stateIds = {
    "RAF", // x_rdata[0]
"pRAF", // x_rdata[1]
"MEK", // x_rdata[2]
"pMEK", // x_rdata[3]
"ERK", // x_rdata[4]
"pERK", // x_rdata[5]
};

std::array<const char*, 12> observableIds = {
    "pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1", // y[0]
"pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2", // y[1]
"pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3", // y[2]
"pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1", // y[3]
"pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2", // y[4]
"pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3", // y[5]
"pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1", // y[6]
"pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2", // y[7]
"pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3", // y[8]
"pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1", // y[9]
"pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2", // y[10]
"pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3", // y[11]
};

std::array<const char*, 7> expressionIds = {
    "k1max", // w[0]
"flux_v1_v_0", // w[1]
"flux_v2_v_1", // w[2]
"flux_v3_v_2", // w[3]
"flux_v4_v_3", // w[4]
"flux_v5_v_4", // w[5]
"flux_v6_v_5", // w[6]
};

std::array<int, 6> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Fiedler_BMCSystBiol2016

} // namespace amici
