#include "wrapfunctions.h"
#include "Laske_PLOSComputBiol2019.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Laske_PLOSComputBiol2019::Model_Laske_PLOSComputBiol2019()
    );
}

} // namespace generic_model

} // namespace amici
