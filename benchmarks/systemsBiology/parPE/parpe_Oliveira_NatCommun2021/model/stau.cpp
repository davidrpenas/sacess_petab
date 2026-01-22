#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "sx.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void stau_Oliveira_NatCommun2021(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 0:
        case 1:
            switch(ip) {
                case 7:
                    stau[0] = -1;
                    break;
            }
            break;
        case 2:
        case 3:
            switch(ip) {
                case 8:
                    stau[0] = -1;
                    break;
            }
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
