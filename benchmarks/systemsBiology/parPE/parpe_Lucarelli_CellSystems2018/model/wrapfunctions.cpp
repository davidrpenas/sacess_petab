#include "wrapfunctions.h"
#include "Lucarelli_CellSystems2018.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Lucarelli_CellSystems2018::Model_Lucarelli_CellSystems2018()
    );
}

} // namespace generic_model

} // namespace amici
