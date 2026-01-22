#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Elowitz_Nature2000 {

// clang-format off

std::array<const char*, 21> parameterNames = {
    "KM", // p[0]
"eff", // p[1]
"eff_GFP", // p[2]
"init_GFP", // p[3]
"init_GFP_mRNA", // p[4]
"init_X_mRNA", // p[5]
"init_X_protein", // p[6]
"init_Y_mRNA", // p[7]
"init_Y_protein", // p[8]
"init_Z_mRNA", // p[9]
"init_Z_protein", // p[10]
"n_Hill", // p[11]
"tau_mRNA", // p[12]
"tau_mRNA_GFP", // p[13]
"tau_prot", // p[14]
"tau_prot_GFP", // p[15]
"tps_active", // p[16]
"tps_repr", // p[17]
"observableParameter1_fluorescence", // p[18]
"observableParameter2_fluorescence", // p[19]
"noiseParameter1_fluorescence", // p[20]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 8> stateNames = {
    "LacI protein", // x_rdata[0]
"TetR protein", // x_rdata[1]
"cI protein", // x_rdata[2]
"GFP", // x_rdata[3]
"LacI mRNA", // x_rdata[4]
"TetR mRNA", // x_rdata[5]
"cI mRNA", // x_rdata[6]
"GFP mRNA", // x_rdata[7]
};

std::array<const char*, 1> observableNames = {
    "", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::log10, // y[0]
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

std::array<const char*, 21> parameterIds = {
    "KM", // p[0]
"eff", // p[1]
"eff_GFP", // p[2]
"init_GFP", // p[3]
"init_GFP_mRNA", // p[4]
"init_X_mRNA", // p[5]
"init_X_protein", // p[6]
"init_Y_mRNA", // p[7]
"init_Y_protein", // p[8]
"init_Z_mRNA", // p[9]
"init_Z_protein", // p[10]
"n_Hill", // p[11]
"tau_mRNA", // p[12]
"tau_mRNA_GFP", // p[13]
"tau_prot", // p[14]
"tau_prot_GFP", // p[15]
"tps_active", // p[16]
"tps_repr", // p[17]
"observableParameter1_fluorescence", // p[18]
"observableParameter2_fluorescence", // p[19]
"noiseParameter1_fluorescence", // p[20]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 8> stateIds = {
    "X_protein", // x_rdata[0]
"Y_protein", // x_rdata[1]
"Z_protein", // x_rdata[2]
"GFP", // x_rdata[3]
"X_mRNA", // x_rdata[4]
"Y_mRNA", // x_rdata[5]
"Z_mRNA", // x_rdata[6]
"GFP_mRNA", // x_rdata[7]
};

std::array<const char*, 1> observableIds = {
    "fluorescence", // y[0]
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

std::array<int, 8> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Elowitz_Nature2000

} // namespace amici
