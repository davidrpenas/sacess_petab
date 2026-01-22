#include "wrapfunctions.h"
#include "Borghans_BiophysChem1997.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Borghans_BiophysChem1997::Model_Borghans_BiophysChem1997()
    );
}

} // namespace generic_model

} // namespace amici
