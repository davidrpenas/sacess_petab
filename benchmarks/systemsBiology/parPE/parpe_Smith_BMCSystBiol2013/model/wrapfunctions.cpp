#include "wrapfunctions.h"
#include "Smith_BMCSystBiol2013.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Smith_BMCSystBiol2013::Model_Smith_BMCSystBiol2013()
    );
}

} // namespace generic_model

} // namespace amici
