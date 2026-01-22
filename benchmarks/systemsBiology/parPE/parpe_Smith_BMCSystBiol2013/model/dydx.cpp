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
#include "dwdx.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void dydx_Smith_BMCSystBiol2013(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[29] = sc_BINS;
    dydx[38] = sc_BINS;
    dydx[47] = sc_BINS;
    dydx[85] = -PTP1B*sc_PTP/std::pow(PTP1B_plus_PTP1B_ox, 2) + sc_PTP/PTP1B_plus_PTP1B_ox;
    dydx[90] = sc_PI3K;
    dydx[181] = sc_GLUT_2B;
    dydx[185] = sc_GLUT_3B_120;
    dydx[186] = sc_GLUT_3B_240;
    dydx[220] = -PTP1B*sc_PTP/std::pow(PTP1B_plus_PTP1B_ox, 2);
    dydx[268] = sc_SOD2;
    dydx[318] = sc_PIRS;
    dydx[431] = sc_FOXO1;
    dydx[440] = sc_FOXO1;
    dydx[449] = sc_FOXO1;
    dydx[458] = sc_FOXO1;
    dydx[467] = sc_FOXO1;
    dydx[476] = sc_FOXO1;
    dydx[485] = sc_FOXO1;
    dydx[494] = sc_FOXO1;
    dydx[503] = sc_FOXO1;
    dydx[512] = sc_FOXO1;
    dydx[521] = sc_FOXO1;
    dydx[530] = sc_FOXO1;
    dydx[539] = sc_FOXO1;
    dydx[548] = sc_FOXO1;
    dydx[557] = sc_FOXO1;
    dydx[566] = sc_FOXO1;
    dydx[575] = sc_FOXO1;
    dydx[584] = sc_FOXO1;
    dydx[593] = sc_FOXO1;
    dydx[602] = sc_FOXO1;
    dydx[611] = sc_FOXO1;
    dydx[620] = sc_FOXO1;
    dydx[629] = sc_FOXO1;
    dydx[638] = sc_FOXO1;
    dydx[647] = sc_FOXO1;
    dydx[656] = sc_FOXO1;
    dydx[665] = sc_FOXO1;
    dydx[674] = sc_FOXO1;
    dydx[683] = sc_FOXO1;
    dydx[692] = sc_FOXO1;
    dydx[701] = sc_FOXO1;
    dydx[710] = sc_FOXO1;
    dydx[719] = sc_FOXO1;
    dydx[728] = sc_FOXO1;
    dydx[737] = sc_FOXO1;
    dydx[746] = sc_FOXO1;
    dydx[755] = sc_FOXO1;
    dydx[764] = sc_FOXO1;
    dydx[773] = sc_FOXO1;
    dydx[782] = sc_FOXO1;
    dydx[791] = sc_FOXO1;
    dydx[800] = sc_FOXO1;
    dydx[809] = sc_FOXO1;
    dydx[818] = sc_FOXO1;
    dydx[827] = sc_FOXO1;
    dydx[836] = sc_FOXO1;
    dydx[845] = sc_FOXO1;
    dydx[854] = sc_FOXO1;
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
