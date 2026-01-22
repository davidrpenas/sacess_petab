#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Perelson_Science1996 {

// clang-format off

std::array<const char*, 3> parameterNames = {
    "c", // p[0]
"delta", // p[1]
"noiseParameter1_task0_model0_perelson1_V", // p[2]
};

std::array<const char*, 3> fixedParameterNames = {
    "NN", // k[0]
"T0", // k[1]
"K0", // k[2]
};

std::array<const char*, 4> stateNames = {
    "Tstar", // x_rdata[0]
"V", // x_rdata[1]
"Vin", // x_rdata[2]
"Vni", // x_rdata[3]
};

std::array<const char*, 1> observableNames = {
    "", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::log10, // y[0]
};

std::array<const char*, 5> expressionNames = {
    "flux_v1", // w[0]
"flux_v2", // w[1]
"flux_v3", // w[2]
"flux_v4", // w[3]
"flux_v5", // w[4]
};

std::array<const char*, 3> parameterIds = {
    "c", // p[0]
"delta", // p[1]
"noiseParameter1_task0_model0_perelson1_V", // p[2]
};

std::array<const char*, 3> fixedParameterIds = {
    "NN", // k[0]
"T0", // k[1]
"K0", // k[2]
};

std::array<const char*, 4> stateIds = {
    "Tstar", // x_rdata[0]
"V", // x_rdata[1]
"Vin", // x_rdata[2]
"Vni", // x_rdata[3]
};

std::array<const char*, 1> observableIds = {
    "task0_model0_perelson1_V", // y[0]
};

std::array<const char*, 5> expressionIds = {
    "flux_v1", // w[0]
"flux_v2", // w[1]
"flux_v3", // w[2]
"flux_v4", // w[3]
"flux_v5", // w[4]
};

std::array<int, 4> stateIdxsSolver = {
    0, 1, 2, 3
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Perelson_Science1996

} // namespace amici
