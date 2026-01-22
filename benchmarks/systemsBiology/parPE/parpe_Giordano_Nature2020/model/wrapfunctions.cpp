#include "wrapfunctions.h"
#include "Giordano_Nature2020.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Giordano_Nature2020::Model_Giordano_Nature2020()
    );
}

} // namespace generic_model

} // namespace amici
