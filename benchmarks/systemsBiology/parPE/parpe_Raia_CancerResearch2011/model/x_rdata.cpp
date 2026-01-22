#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void x_rdata_Raia_CancerResearch2011(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Rec;
    x_rdata[1] = Rec_i;
    x_rdata[2] = IL13_Rec;
    x_rdata[3] = p_IL13_Rec;
    x_rdata[4] = p_IL13_Rec_i;
    x_rdata[5] = JAK2;
    x_rdata[6] = pJAK2;
    x_rdata[7] = STAT5;
    x_rdata[8] = pSTAT5;
    x_rdata[9] = SOCS3mRNA;
    x_rdata[10] = DecoyR;
    x_rdata[11] = IL13_DecoyR;
    x_rdata[12] = SOCS3;
    x_rdata[13] = CD274mRNA;
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
