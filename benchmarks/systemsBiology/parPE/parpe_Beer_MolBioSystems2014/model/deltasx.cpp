#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "w.h"
#include "xdot.h"
#include "xdot_old.h"
#include "sx.h"
#include "stau.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void deltasx_Beer_MolBioSystems2014(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl){
    switch(ie) {
        case 0:
        case 1:
            switch(ip) {
                case 5:
                    deltasx[0] = -dBacdt + xdot_old0;
                    deltasx[1] = -dGludt + xdot_old1;
                    deltasx[2] = -dcGludt + xdot_old2;
                    deltasx[3] = -dInddt + xdot_old3;
                    break;
            }
            break;
    }
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
