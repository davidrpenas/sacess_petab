#include "wrapfunctions.h"
#include "Fiedler_BMCSystBiol2016.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Fiedler_BMCSystBiol2016::Model_Fiedler_BMCSystBiol2016()
    );
}

} // namespace generic_model

} // namespace amici
