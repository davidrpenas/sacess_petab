#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void y_Smith_BMCSystBiol2013(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = IRS1_TyrP_PI3K*sc_PI3K;
    y[1] = cellsurface_GLUT4*sc_GLUT_2B;
    y[2] = InR_bound*sc_BINS;
    y[3] = IRS1_TyrP*sc_PIRS;
    y[4] = PTP1B*sc_PTP/PTP1B_plus_PTP1B_ox;
    y[5] = cellsurface_GLUT4*sc_GLUT_3B_120;
    y[6] = cellsurface_GLUT4*sc_GLUT_3B_240;
    y[7] = cytoplasm_SOD2*sc_SOD2;
    y[8] = Foxo1_all*sc_FOXO1;
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
