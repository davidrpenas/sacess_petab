#include "wrapfunctions.h"
#include "Chen_MSB2009.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Chen_MSB2009::Model_Chen_MSB2009()
    );
}

} // namespace generic_model

} // namespace amici
