#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Bertozzi_PNAS2020 {

// clang-format off

std::array<const char*, 5> parameterNames = {
    "gamma_", // p[0]
"I0_", // p[1]
"R0_", // p[2]
"sd_I_NY", // p[3]
"sd_I_CA", // p[4]
};

std::array<const char*, 1> fixedParameterNames = {
    "N_", // k[0]
};

std::array<const char*, 3> stateNames = {
    "I_", // x_rdata[0]
"R_", // x_rdata[1]
"S_", // x_rdata[2]
};

std::array<const char*, 2> observableNames = {
    "Confirmed infected in NY", // y[0]
"Confirmed infected in CA", // y[1]
};

std::array<const ObservableScaling, 2> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
};

std::array<const char*, 3> expressionNames = {
    "beta_N", // w[0]
"flux_infection", // w[1]
"flux_recovery", // w[2]
};

std::array<const char*, 5> parameterIds = {
    "gamma_", // p[0]
"I0_", // p[1]
"R0_", // p[2]
"sd_I_NY", // p[3]
"sd_I_CA", // p[4]
};

std::array<const char*, 1> fixedParameterIds = {
    "N_", // k[0]
};

std::array<const char*, 3> stateIds = {
    "I_", // x_rdata[0]
"R_", // x_rdata[1]
"S_", // x_rdata[2]
};

std::array<const char*, 2> observableIds = {
    "y_I_NY", // y[0]
"y_I_CA", // y[1]
};

std::array<const char*, 3> expressionIds = {
    "beta_N", // w[0]
"flux_infection", // w[1]
"flux_recovery", // w[2]
};

std::array<int, 3> stateIdxsSolver = {
    0, 1, 2
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Bertozzi_PNAS2020

} // namespace amici
