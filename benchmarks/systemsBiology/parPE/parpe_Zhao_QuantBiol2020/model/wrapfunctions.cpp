#include "wrapfunctions.h"
#include "Zhao_QuantBiol2020.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Zhao_QuantBiol2020::Model_Zhao_QuantBiol2020()
    );
}

} // namespace generic_model

} // namespace amici
