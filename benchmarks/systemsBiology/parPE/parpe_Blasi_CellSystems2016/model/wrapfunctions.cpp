#include "wrapfunctions.h"
#include "Blasi_CellSystems2016.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Blasi_CellSystems2016::Model_Blasi_CellSystems2016()
    );
}

} // namespace generic_model

} // namespace amici
