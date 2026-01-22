#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Crauste_CellSystems2017 {

// clang-format off

std::array<const char*, 33> parameterNames = {
    "delta_EL", // p[0]
"delta_LM", // p[1]
"delta_NE", // p[2]
"mu_EE", // p[3]
"mu_LE", // p[4]
"mu_LL", // p[5]
"mu_N", // p[6]
"mu_P", // p[7]
"mu_PE", // p[8]
"mu_PL", // p[9]
"rho_E", // p[10]
"rho_P", // p[11]
"noiseParameter1_observable_EarlyEffector__266__model1_data1", // p[12]
"noiseParameter1_observable_EarlyEffector__526__model1_data1", // p[13]
"noiseParameter1_observable_EarlyEffector__1540__model1_data1", // p[14]
"noiseParameter1_observable_EarlyEffector__2924__model1_data1", // p[15]
"noiseParameter1_observable_EarlyEffector__3048__model1_data1", // p[16]
"noiseParameter1_observable_EarlyEffector__19722__model1_data1", // p[17]
"noiseParameter1_observable_EarlyEffector__23855__model1_data1", // p[18]
"noiseParameter1_observable_LateEffector__182__model1_data1", // p[19]
"noiseParameter1_observable_LateEffector__1976__model1_data1", // p[20]
"noiseParameter1_observable_LateEffector__2705__model1_data1", // p[21]
"noiseParameter1_observable_LateEffector__4105__model1_data1", // p[22]
"noiseParameter1_observable_LateEffector__6841__model1_data1", // p[23]
"noiseParameter1_observable_LateEffector__7266__model1_data1", // p[24]
"noiseParameter1_observable_LateEffector__23601__model1_data1", // p[25]
"noiseParameter1_observable_Memory__80__model1_data1", // p[26]
"noiseParameter1_observable_Memory__4206__model1_data1", // p[27]
"noiseParameter1_observable_Memory__5745__model1_data1", // p[28]
"noiseParameter1_observable_Memory__7909__model1_data1", // p[29]
"noiseParameter1_observable_Memory__15987__model1_data1", // p[30]
"noiseParameter1_observable_Naive__64__model1_data1", // p[31]
"noiseParameter1_observable_Naive__290__model1_data1", // p[32]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 5> stateNames = {
    "Naive", // x_rdata[0]
"EarlyEffector", // x_rdata[1]
"LateEffector", // x_rdata[2]
"Memory", // x_rdata[3]
"Pathogen", // x_rdata[4]
};

std::array<const char*, 21> observableNames = {
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
"", // y[12]
"", // y[13]
"", // y[14]
"", // y[15]
"", // y[16]
"", // y[17]
"", // y[18]
"", // y[19]
"", // y[20]
};

std::array<const ObservableScaling, 21> observableScalings = {
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
ObservableScaling::lin, // y[12]
ObservableScaling::lin, // y[13]
ObservableScaling::lin, // y[14]
ObservableScaling::lin, // y[15]
ObservableScaling::lin, // y[16]
ObservableScaling::lin, // y[17]
ObservableScaling::lin, // y[18]
ObservableScaling::lin, // y[19]
ObservableScaling::lin, // y[20]
};

std::array<const char*, 12> expressionNames = {
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
};

std::array<const char*, 33> parameterIds = {
    "delta_EL", // p[0]
"delta_LM", // p[1]
"delta_NE", // p[2]
"mu_EE", // p[3]
"mu_LE", // p[4]
"mu_LL", // p[5]
"mu_N", // p[6]
"mu_P", // p[7]
"mu_PE", // p[8]
"mu_PL", // p[9]
"rho_E", // p[10]
"rho_P", // p[11]
"noiseParameter1_observable_EarlyEffector__266__model1_data1", // p[12]
"noiseParameter1_observable_EarlyEffector__526__model1_data1", // p[13]
"noiseParameter1_observable_EarlyEffector__1540__model1_data1", // p[14]
"noiseParameter1_observable_EarlyEffector__2924__model1_data1", // p[15]
"noiseParameter1_observable_EarlyEffector__3048__model1_data1", // p[16]
"noiseParameter1_observable_EarlyEffector__19722__model1_data1", // p[17]
"noiseParameter1_observable_EarlyEffector__23855__model1_data1", // p[18]
"noiseParameter1_observable_LateEffector__182__model1_data1", // p[19]
"noiseParameter1_observable_LateEffector__1976__model1_data1", // p[20]
"noiseParameter1_observable_LateEffector__2705__model1_data1", // p[21]
"noiseParameter1_observable_LateEffector__4105__model1_data1", // p[22]
"noiseParameter1_observable_LateEffector__6841__model1_data1", // p[23]
"noiseParameter1_observable_LateEffector__7266__model1_data1", // p[24]
"noiseParameter1_observable_LateEffector__23601__model1_data1", // p[25]
"noiseParameter1_observable_Memory__80__model1_data1", // p[26]
"noiseParameter1_observable_Memory__4206__model1_data1", // p[27]
"noiseParameter1_observable_Memory__5745__model1_data1", // p[28]
"noiseParameter1_observable_Memory__7909__model1_data1", // p[29]
"noiseParameter1_observable_Memory__15987__model1_data1", // p[30]
"noiseParameter1_observable_Naive__64__model1_data1", // p[31]
"noiseParameter1_observable_Naive__290__model1_data1", // p[32]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 5> stateIds = {
    "Naive", // x_rdata[0]
"EarlyEffector", // x_rdata[1]
"LateEffector", // x_rdata[2]
"Memory", // x_rdata[3]
"Pathogen", // x_rdata[4]
};

std::array<const char*, 21> observableIds = {
    "observable_EarlyEffector__266__model1_data1", // y[0]
"observable_EarlyEffector__526__model1_data1", // y[1]
"observable_EarlyEffector__1540__model1_data1", // y[2]
"observable_EarlyEffector__2924__model1_data1", // y[3]
"observable_EarlyEffector__3048__model1_data1", // y[4]
"observable_EarlyEffector__19722__model1_data1", // y[5]
"observable_EarlyEffector__23855__model1_data1", // y[6]
"observable_LateEffector__182__model1_data1", // y[7]
"observable_LateEffector__1976__model1_data1", // y[8]
"observable_LateEffector__2705__model1_data1", // y[9]
"observable_LateEffector__4105__model1_data1", // y[10]
"observable_LateEffector__6841__model1_data1", // y[11]
"observable_LateEffector__7266__model1_data1", // y[12]
"observable_LateEffector__23601__model1_data1", // y[13]
"observable_Memory__80__model1_data1", // y[14]
"observable_Memory__4206__model1_data1", // y[15]
"observable_Memory__5745__model1_data1", // y[16]
"observable_Memory__7909__model1_data1", // y[17]
"observable_Memory__15987__model1_data1", // y[18]
"observable_Naive__64__model1_data1", // y[19]
"observable_Naive__290__model1_data1", // y[20]
};

std::array<const char*, 12> expressionIds = {
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
};

std::array<int, 5> stateIdxsSolver = {
    0, 1, 2, 3, 4
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Crauste_CellSystems2017

} // namespace amici
