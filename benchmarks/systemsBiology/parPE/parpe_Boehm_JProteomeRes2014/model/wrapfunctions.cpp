#include "wrapfunctions.h"
#include "Boehm_JProteomeRes2014.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Boehm_JProteomeRes2014::Model_Boehm_JProteomeRes2014()
    );
}

} // namespace generic_model

} // namespace amici
