#include "wrapfunctions.h"
#include "Sneyd_PNAS2002.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Sneyd_PNAS2002::Model_Sneyd_PNAS2002()
    );
}

} // namespace generic_model

} // namespace amici
