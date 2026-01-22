#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void xdot_Blasi_CellSystems2016(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dx_0acdt = -flux_re1 - flux_re2 - flux_re3 - flux_re4;  // xdot[0]
    dx_k8dt = -flux_re10 + flux_re2 - flux_re8 - flux_re9;  // xdot[1]
    dx_k5dt = flux_re1 - flux_re5 - flux_re6 - flux_re7;  // xdot[2]
    dx_k12dt = -flux_re11 - flux_re12 - flux_re13 + flux_re3;  // xdot[3]
    dx_k16dt = -flux_re14 - flux_re15 - flux_re16 + flux_re4;  // xdot[4]
    dx_k5k8dt = -flux_re17 - flux_re18 + flux_re5 + flux_re8;  // xdot[5]
    dx_k5k12dt = flux_re11 - flux_re19 - flux_re20 + flux_re6;  // xdot[6]
    dx_k5k16dt = flux_re14 - flux_re21 - flux_re22 + flux_re7;  // xdot[7]
    dx_k8k12dt = flux_re12 - flux_re23 - flux_re24 + flux_re9;  // xdot[8]
    dx_k8k16dt = flux_re10 + flux_re15 - flux_re25 - flux_re26;  // xdot[9]
    dx_k12k16dt = flux_re13 + flux_re16 - flux_re27 - flux_re28;  // xdot[10]
    dx_k5k8k12dt = flux_re17 + flux_re19 + flux_re23 - flux_re29;  // xdot[11]
    dx_k5k8k16dt = flux_re18 + flux_re21 + flux_re25 - flux_re30;  // xdot[12]
    dx_k5k12k16dt = flux_re20 + flux_re22 + flux_re27 - flux_re31;  // xdot[13]
    dx_k8k12k16dt = flux_re24 + flux_re26 + flux_re28 - flux_re32;  // xdot[14]
    dx_4acdt = flux_re29 + flux_re30 + flux_re31 + flux_re32;  // xdot[15]
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
