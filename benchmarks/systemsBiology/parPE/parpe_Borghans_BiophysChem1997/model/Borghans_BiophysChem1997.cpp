#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Borghans_BiophysChem1997 {

// clang-format off

std::array<const char*, 23> parameterNames = {
    "K2", // p[0]
"K_par", // p[1]
"Ka", // p[2]
"Kd", // p[3]
"Kf", // p[4]
"Kp", // p[5]
"Ky", // p[6]
"Kz", // p[7]
"Vd", // p[8]
"Vm2", // p[9]
"Vm3", // p[10]
"Vp", // p[11]
"beta_par", // p[12]
"epsilon_par", // p[13]
"init_A_state", // p[14]
"init_Y_state", // p[15]
"init_Z_state", // p[16]
"n_par", // p[17]
"v0", // p[18]
"v1", // p[19]
"observableParameter1_Ca", // p[20]
"observableParameter2_Ca", // p[21]
"noiseParameter1_Ca", // p[22]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 3> stateNames = {
    "Ca^{2+} in cytosol", // x_rdata[0]
"Ca^{2+} in Ca^{2+}-sensitive pool", // x_rdata[1]
"IP_3 concentration", // x_rdata[2]
};

std::array<const char*, 1> observableNames = {
    "", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::log10, // y[0]
};

std::array<const char*, 8> expressionNames = {
    "flux_v1_v_0", // w[0]
"flux_v2_v_1", // w[1]
"flux_v3_v_2", // w[2]
"flux_v4_v_3", // w[3]
"flux_v5_v_4", // w[4]
"flux_v6_v_5", // w[5]
"flux_v7_v_6", // w[6]
"flux_v8_v_7", // w[7]
};

std::array<const char*, 23> parameterIds = {
    "K2", // p[0]
"K_par", // p[1]
"Ka", // p[2]
"Kd", // p[3]
"Kf", // p[4]
"Kp", // p[5]
"Ky", // p[6]
"Kz", // p[7]
"Vd", // p[8]
"Vm2", // p[9]
"Vm3", // p[10]
"Vp", // p[11]
"beta_par", // p[12]
"epsilon_par", // p[13]
"init_A_state", // p[14]
"init_Y_state", // p[15]
"init_Z_state", // p[16]
"n_par", // p[17]
"v0", // p[18]
"v1", // p[19]
"observableParameter1_Ca", // p[20]
"observableParameter2_Ca", // p[21]
"noiseParameter1_Ca", // p[22]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 3> stateIds = {
    "Z_state", // x_rdata[0]
"Y_state", // x_rdata[1]
"A_state", // x_rdata[2]
};

std::array<const char*, 1> observableIds = {
    "Ca", // y[0]
};

std::array<const char*, 8> expressionIds = {
    "flux_v1_v_0", // w[0]
"flux_v2_v_1", // w[1]
"flux_v3_v_2", // w[2]
"flux_v4_v_3", // w[3]
"flux_v5_v_4", // w[4]
"flux_v6_v_5", // w[5]
"flux_v7_v_6", // w[6]
"flux_v8_v_7", // w[7]
};

std::array<int, 3> stateIdxsSolver = {
    0, 1, 2
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Borghans_BiophysChem1997

} // namespace amici
