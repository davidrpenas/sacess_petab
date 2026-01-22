#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Beer_MolBioSystems2014 {

// clang-format off

std::array<const char*, 9> parameterNames = {
    "Bacmax", // p[0]
"beta", // p[1]
"kdegi", // p[2]
"kdim", // p[3]
"ksyn", // p[4]
"tau", // p[5]
"init_Bac", // p[6]
"noiseParameter1_Bacnorm", // p[7]
"noiseParameter1_IndconcNormRange", // p[8]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 4> stateNames = {
    "Bac", // x_rdata[0]
"Glu", // x_rdata[1]
"cGlu", // x_rdata[2]
"Ind", // x_rdata[3]
};

std::array<const char*, 2> observableNames = {
    "", // y[0]
"", // y[1]
};

std::array<const ObservableScaling, 2> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
};

std::array<const char*, 5> expressionNames = {
    "lag", // w[0]
"flux_reaction1", // w[1]
"flux_reaction2", // w[2]
"flux_reaction3", // w[3]
"flux_reaction4", // w[4]
};

std::array<const char*, 9> parameterIds = {
    "Bacmax", // p[0]
"beta", // p[1]
"kdegi", // p[2]
"kdim", // p[3]
"ksyn", // p[4]
"tau", // p[5]
"init_Bac", // p[6]
"noiseParameter1_Bacnorm", // p[7]
"noiseParameter1_IndconcNormRange", // p[8]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 4> stateIds = {
    "Bac", // x_rdata[0]
"Glu", // x_rdata[1]
"cGlu", // x_rdata[2]
"Ind", // x_rdata[3]
};

std::array<const char*, 2> observableIds = {
    "Bacnorm", // y[0]
"IndconcNormRange", // y[1]
};

std::array<const char*, 5> expressionIds = {
    "lag", // w[0]
"flux_reaction1", // w[1]
"flux_reaction2", // w[2]
"flux_reaction3", // w[3]
"flux_reaction4", // w[4]
};

std::array<int, 4> stateIdxsSolver = {
    0, 1, 2, 3
};

std::array<bool, 2> rootInitialValues = {
    true, true
};

// clang-format on

} // namespace model_Beer_MolBioSystems2014

} // namespace amici
