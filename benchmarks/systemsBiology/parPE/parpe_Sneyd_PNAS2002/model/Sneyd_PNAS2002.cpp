#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Sneyd_PNAS2002 {

// clang-format off

std::array<const char*, 15> parameterNames = {
    "k1", // p[0]
"k2", // p[1]
"k3", // p[2]
"k4", // p[3]
"k_1", // p[4]
"k_2", // p[5]
"k_3", // p[6]
"k_4", // p[7]
"l2", // p[8]
"l4", // p[9]
"l6", // p[10]
"l_2", // p[11]
"l_4", // p[12]
"l_6", // p[13]
"noiseParameter1_open_probability", // p[14]
};

std::array<const char*, 2> fixedParameterNames = {
    "Ca", // k[0]
"IP3", // k[1]
};

std::array<const char*, 6> stateNames = {
    "O", // x_rdata[0]
"R", // x_rdata[1]
"I_1", // x_rdata[2]
"S", // x_rdata[3]
"A", // x_rdata[4]
"I_2", // x_rdata[5]
};

std::array<const char*, 1> observableNames = {
    "", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 13> expressionNames = {
    "L1", // w[0]
"L3", // w[1]
"L5", // w[2]
"flux_v_0", // w[3]
"flux_v_1", // w[4]
"flux_v_2", // w[5]
"flux_v_3", // w[6]
"flux_v_4", // w[7]
"flux_v_5", // w[8]
"flux_v_6", // w[9]
"flux_v_7", // w[10]
"flux_v_8", // w[11]
"flux_v_9", // w[12]
};

std::array<const char*, 15> parameterIds = {
    "k1", // p[0]
"k2", // p[1]
"k3", // p[2]
"k4", // p[3]
"k_1", // p[4]
"k_2", // p[5]
"k_3", // p[6]
"k_4", // p[7]
"l2", // p[8]
"l4", // p[9]
"l6", // p[10]
"l_2", // p[11]
"l_4", // p[12]
"l_6", // p[13]
"noiseParameter1_open_probability", // p[14]
};

std::array<const char*, 2> fixedParameterIds = {
    "Ca", // k[0]
"IP3", // k[1]
};

std::array<const char*, 6> stateIds = {
    "IPR_O", // x_rdata[0]
"IPR_R", // x_rdata[1]
"IPR_I1", // x_rdata[2]
"IPR_S", // x_rdata[3]
"IPR_A", // x_rdata[4]
"IPR_I2", // x_rdata[5]
};

std::array<const char*, 1> observableIds = {
    "open_probability", // y[0]
};

std::array<const char*, 13> expressionIds = {
    "L1", // w[0]
"L3", // w[1]
"L5", // w[2]
"flux_v_0", // w[3]
"flux_v_1", // w[4]
"flux_v_2", // w[5]
"flux_v_3", // w[6]
"flux_v_4", // w[7]
"flux_v_5", // w[8]
"flux_v_6", // w[9]
"flux_v_7", // w[10]
"flux_v_8", // w[11]
"flux_v_9", // w[12]
};

std::array<int, 6> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Sneyd_PNAS2002

} // namespace amici
