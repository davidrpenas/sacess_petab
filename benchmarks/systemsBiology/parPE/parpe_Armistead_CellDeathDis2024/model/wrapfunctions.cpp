#include "wrapfunctions.h"
#include "Armistead_CellDeathDis2024.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Armistead_CellDeathDis2024::Model_Armistead_CellDeathDis2024()
    );
}

} // namespace generic_model

} // namespace amici
