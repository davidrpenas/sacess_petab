#include "wrapfunctions.h"
#include "Beer_MolBioSystems2014.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Beer_MolBioSystems2014::Model_Beer_MolBioSystems2014()
    );
}

} // namespace generic_model

} // namespace amici
