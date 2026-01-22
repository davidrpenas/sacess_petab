#include <amici/defines.h>
#include <array>

namespace amici {

namespace model_Rahman_MBS2016 {

// clang-format off

std::array<const char*, 10> parameterNames = {
    "infected_normal_transmission_rate_relative", // p[0]
"infected_moderate_transmission_rate", // p[1]
"infected_weak_transmission_rate_relative", // p[2]
"infected_weak_treatment_rate", // p[3]
"infected_normal_worsen_rate", // p[4]
"infected_moderate_worsen_rate", // p[5]
"treated_moderate_improve_rate", // p[6]
"treated_weak_improve_rate", // p[7]
"behavioural_change_rate", // p[8]
"noiseParameter1_observable_prevalence", // p[9]
};

std::array<const char*, 10> fixedParameterNames = {
    "recruitment_rate", // k[0]
"infected_normal_treatment_rate", // k[1]
"infected_moderate_treatment_rate", // k[2]
"susceptible_death_rate", // k[3]
"infected_normal_death_rate", // k[4]
"infected_moderate_death_rate", // k[5]
"infected_weak_death_rate", // k[6]
"treated_normal_death_rate", // k[7]
"treated_moderate_death_rate", // k[8]
"treated_weak_death_rate", // k[9]
};

std::array<const char*, 7> stateNames = {
    "susceptible", // x_rdata[0]
"infected_normal", // x_rdata[1]
"infected_moderate", // x_rdata[2]
"infected_weak", // x_rdata[3]
"treated_normal", // x_rdata[4]
"treated_moderate", // x_rdata[5]
"treated_weak", // x_rdata[6]
};

std::array<const char*, 1> observableNames = {
    "prevalence", // y[0]
};

std::array<const ObservableScaling, 1> observableScalings = {
    ObservableScaling::lin, // y[0]
};

std::array<const char*, 21> expressionNames = {
    "infected_normal_transmission_rate", // w[0]
"infected_weak_transmission_rate", // w[1]
"prevalence", // w[2]
"treated_transmission_rate", // w[3]
"force_of_infection", // w[4]
"flux_v1", // w[5]
"flux_v2", // w[6]
"flux_v3", // w[7]
"flux_v4", // w[8]
"flux_v5", // w[9]
"flux_v6", // w[10]
"flux_v7", // w[11]
"flux_v8", // w[12]
"flux_v9", // w[13]
"flux_v10", // w[14]
"flux_v11", // w[15]
"flux_v12", // w[16]
"flux_v13", // w[17]
"flux_v14", // w[18]
"flux_v15", // w[19]
"flux_v16", // w[20]
};

std::array<const char*, 10> parameterIds = {
    "infected_normal_transmission_rate_relative", // p[0]
"infected_moderate_transmission_rate", // p[1]
"infected_weak_transmission_rate_relative", // p[2]
"infected_weak_treatment_rate", // p[3]
"infected_normal_worsen_rate", // p[4]
"infected_moderate_worsen_rate", // p[5]
"treated_moderate_improve_rate", // p[6]
"treated_weak_improve_rate", // p[7]
"behavioural_change_rate", // p[8]
"noiseParameter1_observable_prevalence", // p[9]
};

std::array<const char*, 10> fixedParameterIds = {
    "recruitment_rate", // k[0]
"infected_normal_treatment_rate", // k[1]
"infected_moderate_treatment_rate", // k[2]
"susceptible_death_rate", // k[3]
"infected_normal_death_rate", // k[4]
"infected_moderate_death_rate", // k[5]
"infected_weak_death_rate", // k[6]
"treated_normal_death_rate", // k[7]
"treated_moderate_death_rate", // k[8]
"treated_weak_death_rate", // k[9]
};

std::array<const char*, 7> stateIds = {
    "susceptible", // x_rdata[0]
"infected_normal", // x_rdata[1]
"infected_moderate", // x_rdata[2]
"infected_weak", // x_rdata[3]
"treated_normal", // x_rdata[4]
"treated_moderate", // x_rdata[5]
"treated_weak", // x_rdata[6]
};

std::array<const char*, 1> observableIds = {
    "observable_prevalence", // y[0]
};

std::array<const char*, 21> expressionIds = {
    "infected_normal_transmission_rate", // w[0]
"infected_weak_transmission_rate", // w[1]
"prevalence", // w[2]
"treated_transmission_rate", // w[3]
"force_of_infection", // w[4]
"flux_v1", // w[5]
"flux_v2", // w[6]
"flux_v3", // w[7]
"flux_v4", // w[8]
"flux_v5", // w[9]
"flux_v6", // w[10]
"flux_v7", // w[11]
"flux_v8", // w[12]
"flux_v9", // w[13]
"flux_v10", // w[14]
"flux_v11", // w[15]
"flux_v12", // w[16]
"flux_v13", // w[17]
"flux_v14", // w[18]
"flux_v15", // w[19]
"flux_v16", // w[20]
};

std::array<int, 7> stateIdxsSolver = {
    0, 1, 2, 3, 4, 5, 6
};

std::array<bool, 0> rootInitialValues = {
    
};

// clang-format on

} // namespace model_Rahman_MBS2016

} // namespace amici
