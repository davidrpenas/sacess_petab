#include "wrapfunctions.h"
#include "Fujita_SciSignal2010.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Fujita_SciSignal2010::Model_Fujita_SciSignal2010()
    );
}

} // namespace generic_model

} // namespace amici
