#include "wrapfunctions.h"
#include "Weber_BMC2015.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Weber_BMC2015::Model_Weber_BMC2015()
    );
}

} // namespace generic_model

} // namespace amici
