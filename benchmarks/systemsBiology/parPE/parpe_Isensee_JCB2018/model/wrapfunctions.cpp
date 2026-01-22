#include "wrapfunctions.h"
#include "Isensee_JCB2018.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Isensee_JCB2018::Model_Isensee_JCB2018()
    );
}

} // namespace generic_model

} // namespace amici
