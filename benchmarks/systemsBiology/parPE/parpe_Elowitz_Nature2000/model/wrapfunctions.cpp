#include "wrapfunctions.h"
#include "Elowitz_Nature2000.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Elowitz_Nature2000::Model_Elowitz_Nature2000()
    );
}

} // namespace generic_model

} // namespace amici
