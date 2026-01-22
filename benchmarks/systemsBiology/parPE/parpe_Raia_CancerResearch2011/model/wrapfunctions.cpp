#include "wrapfunctions.h"
#include "Raia_CancerResearch2011.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Raia_CancerResearch2011::Model_Raia_CancerResearch2011()
    );
}

} // namespace generic_model

} // namespace amici
