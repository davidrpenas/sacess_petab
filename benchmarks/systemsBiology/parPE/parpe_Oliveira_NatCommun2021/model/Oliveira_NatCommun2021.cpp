#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Oliveira_NatCommun2021 {

// clang-format off

std::array<const char*, 14> parameterNames = {
    "delta_", // p[0]
"h_hosp_rate", // p[1]
"gamma_h", // p[2]
"gamma_u", // p[3]
"beta_0", // p[4]
"beta_1", // p[5]
"beta_2", // p[6]
"t_1", // p[7]
"t_2", // p[8]
"exposed_init_concentration", // p[9]
"asymptomatic_init_concentration", // p[10]
"symptomatic_init_concentration", // p[11]
"noiseParameter1_cumulative_deaths", // p[12]
"noiseParameter1_cumulative_cases", // p[13]
};

std::array<const char*, 11> fixedParameterNames = {
    "population", // k[0]
"p_symp_rate", // k[1]
"kappa", // k[2]
"gamma_a", // k[3]
"gamma_s", // k[4]
"xi", // k[5]
"mu_h", // k[6]
"mu_u", // k[7]
"omega_h", // k[8]
"omega_u", // k[9]
"beta_2_multiplier", // k[10]
};

std::array<const char*, 9> stateNames = {
    "Exposed", // x_rdata[0]
"Susceptible", // x_rdata[1]
"Asymptomatic", // x_rdata[2]
"Symptomatic", // x_rdata[3]
"Recovered", // x_rdata[4]
"Hospital", // x_rdata[5]
"ICU", // x_rdata[6]
"Deaths", // x_rdata[7]
"Cumulative_cases", // x_rdata[8]
};

std::array<const char*, 2> observableNames = {
    "cumulative_deaths", // y[0]
"cumulative_cases", // y[1]
};

std::array<const ObservableScaling, 2> observableScalings = {
    ObservableScaling::lin, // y[0]
ObservableScaling::lin, // y[1]
};

std::array<const char*, 22> expressionNames = {
    "beta", // w[0]
"exposed_to_asymptomatic", // w[1]
"exposed_to_symptomatic", // w[2]
"hospitalized_to_deaths", // w[3]
"hospitalized_to_icu", // w[4]
"hospitalized_to_recovered", // w[5]
"icu_to_deaths", // w[6]
"symptomatic_to_hospitalized", // w[7]
"symptomatic_to_icu", // w[8]
"symptomatic_to_recovered", // w[9]
"flux_Susceptible_to_Exposed", // w[10]
"flux_Exposed_to_Asymptomatic", // w[11]
"flux_Exposed_to_Symptomatic", // w[12]
"flux_Asymptomatic_to_Recovered", // w[13]
"flux_Symptomatic_to_Recovered", // w[14]
"flux_Symptomatic_to_Hospitalized", // w[15]
"flux_Hospitalized_to_Recovered", // w[16]
"flux_Symptomatic_to_ICU", // w[17]
"flux_Hospitalized_to_ICU", // w[18]
"flux_ICU_to_Hospitalized", // w[19]
"flux_Hospitalized_to_Deaths", // w[20]
"flux_ICU_to_Deaths", // w[21]
};

std::array<const char*, 14> parameterIds = {
    "delta_", // p[0]
"h_hosp_rate", // p[1]
"gamma_h", // p[2]
"gamma_u", // p[3]
"beta_0", // p[4]
"beta_1", // p[5]
"beta_2", // p[6]
"t_1", // p[7]
"t_2", // p[8]
"exposed_init_concentration", // p[9]
"asymptomatic_init_concentration", // p[10]
"symptomatic_init_concentration", // p[11]
"noiseParameter1_cumulative_deaths", // p[12]
"noiseParameter1_cumulative_cases", // p[13]
};

std::array<const char*, 11> fixedParameterIds = {
    "population", // k[0]
"p_symp_rate", // k[1]
"kappa", // k[2]
"gamma_a", // k[3]
"gamma_s", // k[4]
"xi", // k[5]
"mu_h", // k[6]
"mu_u", // k[7]
"omega_h", // k[8]
"omega_u", // k[9]
"beta_2_multiplier", // k[10]
};

std::array<const char*, 9> stateIds = {
    "Exposed", // x_rdata[0]
"Susceptible", // x_rdata[1]
"Asymptomatic", // x_rdata[2]
"Symptomatic", // x_rdata[3]
"Recovered", // x_rdata[4]
"Hospital", // x_rdata[5]
"ICU", // x_rdata[6]
"Deaths", // x_rdata[7]
"Cumulative_cases", // x_rdata[8]
};

std::array<const char*, 2> observableIds = {
    "cumulative_deaths", // y[0]
"cumulative_cases", // y[1]
};

std::array<const char*, 22> expressionIds = {
    "beta", // w[0]
"exposed_to_asymptomatic", // w[1]
"exposed_to_symptomatic", // w[2]
"hospitalized_to_deaths", // w[3]
"hospitalized_to_icu", // w[4]
"hospitalized_to_recovered", // w[5]
"icu_to_deaths", // w[6]
"symptomatic_to_hospitalized", // w[7]
"symptomatic_to_icu", // w[8]
"symptomatic_to_recovered", // w[9]
"flux_Susceptible_to_Exposed", // w[10]
"flux_Exposed_to_Asymptomatic", // w[11]
"flux_Exposed_to_Symptomatic", // w[12]
"flux_Asymptomatic_to_Recovered", // w[13]
"flux_Symptomatic_to_Recovered", // w[14]
"flux_Symptomatic_to_Hospitalized", // w[15]
"flux_Hospitalized_to_Recovered", // w[16]
"flux_Symptomatic_to_ICU", // w[17]
"flux_Hospitalized_to_ICU", // w[18]
"flux_ICU_to_Hospitalized", // w[19]
"flux_Hospitalized_to_Deaths", // w[20]
"flux_ICU_to_Deaths", // w[21]
};

std::array<int, 9> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6, 7, 8
};

std::array<bool, 4> rootInitialValues = {
    true, true, true, true
};

// clang-format on

} // namespace model_Oliveira_NatCommun2021

} // namespace amici
