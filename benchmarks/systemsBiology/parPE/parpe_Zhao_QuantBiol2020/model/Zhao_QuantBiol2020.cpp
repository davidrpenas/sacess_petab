#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Zhao_QuantBiol2020 {

// clang-format off

std::array<const char*, 25> parameterNames = {
    "R_Stage_I_Wuhan", // p[0]
"gamma_1_Stage_I_Wuhan", // p[1]
"gamma_2_Stage_I_Wuhan", // p[2]
"R_Stage_II_Wuhan", // p[3]
"gamma_1_Stage_II_Wuhan", // p[4]
"gamma_2_Stage_II_Wuhan", // p[5]
"R_Stage_III_Wuhan", // p[6]
"gamma_1_Stage_III_Wuhan", // p[7]
"gamma_2_Stage_III_Wuhan", // p[8]
"R_Stage_I_Hubei", // p[9]
"gamma_1_Stage_I_Hubei", // p[10]
"gamma_2_Stage_I_Hubei", // p[11]
"R_Stage_II_Hubei", // p[12]
"gamma_1_Stage_II_Hubei", // p[13]
"gamma_2_Stage_II_Hubei", // p[14]
"R_Stage_I_China", // p[15]
"gamma_1_Stage_I_China", // p[16]
"gamma_2_Stage_I_China", // p[17]
"R_Stage_II_China", // p[18]
"gamma_1_Stage_II_China", // p[19]
"gamma_2_Stage_II_China", // p[20]
"noiseParameter1_observable_confirmed_infected", // p[21]
"noiseParameter1_observable_susceptible", // p[22]
"noiseParameter1_observable_quarantined_infected", // p[23]
"noiseParameter1_observable_cumulative_infected", // p[24]
};

std::array<const char*, 10> fixedParameterNames = {
    "Trigger_Stage_I", // k[0]
"Trigger_Stage_II", // k[1]
"Trigger_Stage_III", // k[2]
"Trigger_Wuhan", // k[3]
"Trigger_Hubei", // k[4]
"Trigger_China", // k[5]
"sigma", // k[6]
"Total_Pop_Wuhan", // k[7]
"Total_Pop_Hubei", // k[8]
"Total_Pop_China", // k[9]
};

std::array<const char*, 4> stateNames = {
    "Susceptible", // x_rdata[0]
"Unquarantined_Infected", // x_rdata[1]
"Quarantined_Infected", // x_rdata[2]
"Confirmed_Infected", // x_rdata[3]
};

std::array<const char*, 5> observableNames = {
    "confirmed_infected", // y[0]
"susceptible", // y[1]
"unquarantined_infected", // y[2]
"quarantined_infected", // y[3]
"cumulative_infected", // y[4]
};

std::array<const ObservableScaling, 5> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
ObservableScaling::lin, // y[2]
ObservableScaling::lin, // y[3]
ObservableScaling::lin, // y[4]
};

std::array<const char*, 10> expressionNames = {
    "Cumulative_Infected", // w[0]
"R", // w[1]
"Total_Pop", // w[2]
"gamma_1", // w[3]
"gamma_2", // w[4]
"alpha", // w[5]
"beta", // w[6]
"flux_Susceptible_to_Unquarantined", // w[7]
"flux_Unquarantined_to_Quarantined", // w[8]
"flux_Quarantined_to_Confirmed", // w[9]
};

std::array<const char*, 25> parameterIds = {
    "R_Stage_I_Wuhan", // p[0]
"gamma_1_Stage_I_Wuhan", // p[1]
"gamma_2_Stage_I_Wuhan", // p[2]
"R_Stage_II_Wuhan", // p[3]
"gamma_1_Stage_II_Wuhan", // p[4]
"gamma_2_Stage_II_Wuhan", // p[5]
"R_Stage_III_Wuhan", // p[6]
"gamma_1_Stage_III_Wuhan", // p[7]
"gamma_2_Stage_III_Wuhan", // p[8]
"R_Stage_I_Hubei", // p[9]
"gamma_1_Stage_I_Hubei", // p[10]
"gamma_2_Stage_I_Hubei", // p[11]
"R_Stage_II_Hubei", // p[12]
"gamma_1_Stage_II_Hubei", // p[13]
"gamma_2_Stage_II_Hubei", // p[14]
"R_Stage_I_China", // p[15]
"gamma_1_Stage_I_China", // p[16]
"gamma_2_Stage_I_China", // p[17]
"R_Stage_II_China", // p[18]
"gamma_1_Stage_II_China", // p[19]
"gamma_2_Stage_II_China", // p[20]
"noiseParameter1_observable_confirmed_infected", // p[21]
"noiseParameter1_observable_susceptible", // p[22]
"noiseParameter1_observable_quarantined_infected", // p[23]
"noiseParameter1_observable_cumulative_infected", // p[24]
};

std::array<const char*, 10> fixedParameterIds = {
    "Trigger_Stage_I", // k[0]
"Trigger_Stage_II", // k[1]
"Trigger_Stage_III", // k[2]
"Trigger_Wuhan", // k[3]
"Trigger_Hubei", // k[4]
"Trigger_China", // k[5]
"sigma", // k[6]
"Total_Pop_Wuhan", // k[7]
"Total_Pop_Hubei", // k[8]
"Total_Pop_China", // k[9]
};

std::array<const char*, 4> stateIds = {
    "Susceptible", // x_rdata[0]
"Unquarantined_Infected", // x_rdata[1]
"Quarantined_Infected", // x_rdata[2]
"Confirmed_Infected", // x_rdata[3]
};

std::array<const char*, 5> observableIds = {
    "observable_confirmed_infected", // y[0]
"observable_susceptible", // y[1]
"observable_unquarantined_infected", // y[2]
"observable_quarantined_infected", // y[3]
"observable_cumulative_infected", // y[4]
};

std::array<const char*, 10> expressionIds = {
    "Cumulative_Infected", // w[0]
"R", // w[1]
"Total_Pop", // w[2]
"gamma_1", // w[3]
"gamma_2", // w[4]
"alpha", // w[5]
"beta", // w[6]
"flux_Susceptible_to_Unquarantined", // w[7]
"flux_Unquarantined_to_Quarantined", // w[8]
"flux_Quarantined_to_Confirmed", // w[9]
};

std::array<int, 4> stateIdxsSolver = {
    0, 1, 2, 3
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Zhao_QuantBiol2020

} // namespace amici
