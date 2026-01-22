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
namespace model_Isensee_JCB2018 {

void root_Isensee_JCB2018(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = -Fsk_time + t;
    root[1] = Fsk_time - t;
    root[2] = -H89_time + t;
    root[3] = H89_time - t;
    root[4] = -IBMX_time + t;
    root[5] = IBMX_time - t;
    root[6] = -Rp8_Br_cAMPS_pAB_incubation_time + t;
    root[7] = Rp8_Br_cAMPS_pAB_incubation_time - t;
    root[8] = -Rp8_pCPT_cAMPS_pAB_incubation_time + t;
    root[9] = Rp8_pCPT_cAMPS_pAB_incubation_time - t;
    root[10] = -Rp_cAMPS_pAB_incubation_time + t;
    root[11] = Rp_cAMPS_pAB_incubation_time - t;
    root[12] = -Sp8_Br_cAMPS_AM_time + t;
    root[13] = Sp8_Br_cAMPS_AM_time - t;
    root[14] = -fourABnOH_incubation_time + t;
    root[15] = fourABnOH_incubation_time - t;
}

} // namespace model_Isensee_JCB2018
} // namespace amici
