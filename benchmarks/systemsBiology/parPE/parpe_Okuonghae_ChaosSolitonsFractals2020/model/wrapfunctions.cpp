#include "wrapfunctions.h"
#include "Okuonghae_ChaosSolitonsFractals2020.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Okuonghae_ChaosSolitonsFractals2020::Model_Okuonghae_ChaosSolitonsFractals2020()
    );
}

} // namespace generic_model

} // namespace amici
