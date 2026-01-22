#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Okuonghae_ChaosSolitonsFractals2020 {

// clang-format off

std::array<const char*, 16> parameterNames = {
    "transmission_rate_effective", // p[0]
"gamma_a", // p[1]
"gamma_0", // p[2]
"gamma_i", // p[3]
"psi", // p[4]
"theta", // p[5]
"d_0", // p[6]
"d_D", // p[7]
"sigma", // p[8]
"nu", // p[9]
"alpha", // p[10]
"exposed_start", // p[11]
"asymptomatic_start", // p[12]
"symptomatic_start", // p[13]
"noiseParameter1_observable_Cumulative", // p[14]
"noiseParameter1_observable_ActiveCases", // p[15]
};

std::array<const char*, 2> fixedParameterNames = {
    "eps", // k[0]
"delta", // k[1]
};

std::array<const char*, 8> stateNames = {
    "susceptible", // x_rdata[0]
"exposed", // x_rdata[1]
"asymptomatic", // x_rdata[2]
"symptomatic", // x_rdata[3]
"detected", // x_rdata[4]
"recovered", // x_rdata[5]
"deceased", // x_rdata[6]
"detected_cumulative", // x_rdata[7]
};

std::array<const char*, 2> observableNames = {
    "detected_cumulative", // y[0]
"ActiveCases", // y[1]
};

std::array<const ObservableScaling, 2> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
};

std::array<const char*, 12> expressionNames = {
    "total_pop", // w[0]
"lam", // w[1]
"flux_susceptible_to_exposed", // w[2]
"flux_exposed_to_symptomatic", // w[3]
"flux_exposed_to_asymptopmatic", // w[4]
"flux_symptomatic_to_detected", // w[5]
"flux_asymptomatic_to_detected", // w[6]
"flux_symptomatic_to_recovered", // w[7]
"flux_asymptomatic_to_recovered", // w[8]
"flux_detected_to_recovered", // w[9]
"flux_symptomatic_to_deceased", // w[10]
"flux_detected_to_deceased", // w[11]
};

std::array<const char*, 16> parameterIds = {
    "transmission_rate_effective", // p[0]
"gamma_a", // p[1]
"gamma_0", // p[2]
"gamma_i", // p[3]
"psi", // p[4]
"theta", // p[5]
"d_0", // p[6]
"d_D", // p[7]
"sigma", // p[8]
"nu", // p[9]
"alpha", // p[10]
"exposed_start", // p[11]
"asymptomatic_start", // p[12]
"symptomatic_start", // p[13]
"noiseParameter1_observable_Cumulative", // p[14]
"noiseParameter1_observable_ActiveCases", // p[15]
};

std::array<const char*, 2> fixedParameterIds = {
    "eps", // k[0]
"delta", // k[1]
};

std::array<const char*, 8> stateIds = {
    "susceptible", // x_rdata[0]
"exposed", // x_rdata[1]
"asymptomatic", // x_rdata[2]
"symptomatic", // x_rdata[3]
"detected", // x_rdata[4]
"recovered", // x_rdata[5]
"deceased", // x_rdata[6]
"detected_cumulative", // x_rdata[7]
};

std::array<const char*, 2> observableIds = {
    "observable_Cumulative", // y[0]
"observable_ActiveCases", // y[1]
};

std::array<const char*, 12> expressionIds = {
    "total_pop", // w[0]
"lam", // w[1]
"flux_susceptible_to_exposed", // w[2]
"flux_exposed_to_symptomatic", // w[3]
"flux_exposed_to_asymptopmatic", // w[4]
"flux_symptomatic_to_detected", // w[5]
"flux_asymptomatic_to_detected", // w[6]
"flux_symptomatic_to_recovered", // w[7]
"flux_asymptomatic_to_recovered", // w[8]
"flux_detected_to_recovered", // w[9]
"flux_symptomatic_to_deceased", // w[10]
"flux_detected_to_deceased", // w[11]
};

std::array<int, 8> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Okuonghae_ChaosSolitonsFractals2020

} // namespace amici
