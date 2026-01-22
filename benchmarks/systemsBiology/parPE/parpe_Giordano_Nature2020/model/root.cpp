#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void root_Giordano_Nature2020(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = initialTimeManual - t + 4;
    root[1] = -initialTimeManual + t - 4;
    root[2] = initialTimeManual - t + 50;
    root[3] = -initialTimeManual + t - 50;
    root[4] = initialTimeManual - t + 22;
    root[5] = -initialTimeManual + t - 22;
    root[6] = initialTimeManual - t + 28;
    root[7] = -initialTimeManual + t - 28;
    root[8] = initialTimeManual - t + 12;
    root[9] = -initialTimeManual + t - 12;
    root[10] = initialTimeManual - t + 38;
    root[11] = -initialTimeManual + t - 38;
}

} // namespace model_Giordano_Nature2020
} // namespace amici
