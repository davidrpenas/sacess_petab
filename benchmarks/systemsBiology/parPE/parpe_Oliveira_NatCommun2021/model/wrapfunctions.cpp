#include "wrapfunctions.h"
#include "Oliveira_NatCommun2021.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Oliveira_NatCommun2021::Model_Oliveira_NatCommun2021()
    );
}

} // namespace generic_model

} // namespace amici
