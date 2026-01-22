#include "wrapfunctions.h"
#include "Perelson_Science1996.h"
#include "amici/model.h"

namespace amici {
namespace generic_model {

std::unique_ptr<amici::Model> getModel() {
    return std::unique_ptr<amici::Model>(
        new amici::model_Perelson_Science1996::Model_Perelson_Science1996()
    );
}

} // namespace generic_model

} // namespace amici
