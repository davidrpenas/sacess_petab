#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Armistead_CellDeathDis2024 {

// clang-format off

std::array<const char*, 18> parameterNames = {
    "Basal Sphinganine production rate", // p[0]
"Sphinganine to Ceramide conversion rate", // p[1]
"Ceramide to Sphingosine conversion rate", // p[2]
"Sphingosine to Ceramide conversion rate", // p[3]
"Sphingosine to S1P conversion rate", // p[4]
"S1P to Sphingosine conversion rate", // p[5]
"S1P degradation rate", // p[6]
"Ceramide degradation rate", // p[7]
"feedback factor", // p[8]
"Rate of Hai1a inhibition", // p[9]
"noiseParameter1_Sphinga_obs", // p[10]
"noiseParameter1_Cer_obs", // p[11]
"noiseParameter1_Sphingo_obs", // p[12]
"noiseParameter1_S1P_obs", // p[13]
"initial_Cer_sim", // p[14]
"initial_Sphinga_sim", // p[15]
"initial_Sphingo_sim", // p[16]
"initial_S1P_sim", // p[17]
};

std::array<const char*, 1> fixedParameterNames = {
    "Boolean Parameter if Hai1a is active or not", // k[0]
};

std::array<const char*, 4> stateNames = {
    "Ceramide", // x_rdata[0]
"Sphingosine", // x_rdata[1]
"Sphingosine-1-Phosphate", // x_rdata[2]
"Sphinganine", // x_rdata[3]
};

std::array<const char*, 4> observableNames = {
    "Sphinganine concentration [area ratio/Protein] - indication of Cer de novo production", // y[0]
"Ceramide concentration [area ratio/Protein]", // y[1]
"Sphingosine concentration [area ratio/Protein]", // y[2]
"Sphingosine-1-phosphate indirect measurement", // y[3]
};

std::array<const ObservableScaling, 4> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
ObservableScaling::lin, // y[3]
};

std::array<const char*, 8> expressionNames = {
    "flux_R1", // w[0]
"flux_R2", // w[1]
"flux_R3", // w[2]
"flux_R4", // w[3]
"flux_R5", // w[4]
"flux_R6", // w[5]
"flux_R7", // w[6]
"flux_R8", // w[7]
};

std::array<const char*, 18> parameterIds = {
    "k00", // p[0]
"k0", // p[1]
"k1", // p[2]
"k2", // p[3]
"k3", // p[4]
"k4", // p[5]
"k5", // p[6]
"k_d", // p[7]
"alpha_cer", // p[8]
"alpha_hai1a", // p[9]
"noiseParameter1_Sphinga_obs", // p[10]
"noiseParameter1_Cer_obs", // p[11]
"noiseParameter1_Sphingo_obs", // p[12]
"noiseParameter1_S1P_obs", // p[13]
"initial_Cer_sim", // p[14]
"initial_Sphinga_sim", // p[15]
"initial_Sphingo_sim", // p[16]
"initial_S1P_sim", // p[17]
};

std::array<const char*, 1> fixedParameterIds = {
    "S_on", // k[0]
};

std::array<const char*, 4> stateIds = {
    "Cer", // x_rdata[0]
"Sphingo", // x_rdata[1]
"S1P", // x_rdata[2]
"Sphinga", // x_rdata[3]
};

std::array<const char*, 4> observableIds = {
    "Sphinga_obs", // y[0]
"Cer_obs", // y[1]
"Sphingo_obs", // y[2]
"S1P_obs", // y[3]
};

std::array<const char*, 8> expressionIds = {
    "flux_R1", // w[0]
"flux_R2", // w[1]
"flux_R3", // w[2]
"flux_R4", // w[3]
"flux_R5", // w[4]
"flux_R6", // w[5]
"flux_R7", // w[6]
"flux_R8", // w[7]
};

std::array<int, 4> stateIdxsSolver = {
    0, 1, 2, 3
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Armistead_CellDeathDis2024

} // namespace amici
