#include "wrapfunctions.h"
#include "Schwen_PONE2014.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Schwen_PONE2014::Model_Schwen_PONE2014()
    );
}

} // namespace generic_model

} // namespace amici
