#include "wrapfunctions.h"
#include "Crauste_CellSystems2017.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Crauste_CellSystems2017::Model_Crauste_CellSystems2017()
    );
}

} // namespace generic_model

} // namespace amici
