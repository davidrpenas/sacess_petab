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
namespace model_Brannmark_JBC2010 {

void root_Brannmark_JBC2010(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = -insulin_time_1 + t;
    root[1] = insulin_time_1 - t;
    root[2] = -insulin_time_2 + t;
    root[3] = insulin_time_2 - t;
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
