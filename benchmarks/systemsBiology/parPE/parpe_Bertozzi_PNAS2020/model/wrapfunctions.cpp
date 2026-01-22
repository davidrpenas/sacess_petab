#include "wrapfunctions.h"
#include "Bertozzi_PNAS2020.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Bertozzi_PNAS2020::Model_Bertozzi_PNAS2020()
    );
}

} // namespace generic_model

} // namespace amici
