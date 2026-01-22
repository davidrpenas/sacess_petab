#include "wrapfunctions.h"
#include "Zheng_PNAS2012.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Zheng_PNAS2012::Model_Zheng_PNAS2012()
    );
}

} // namespace generic_model

} // namespace amici
