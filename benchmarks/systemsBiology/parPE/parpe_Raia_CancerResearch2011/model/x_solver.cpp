#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void x_solver_Raia_CancerResearch2011(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Rec;
    x_solver[1] = Rec_i;
    x_solver[2] = IL13_Rec;
    x_solver[3] = p_IL13_Rec;
    x_solver[4] = p_IL13_Rec_i;
    x_solver[5] = JAK2;
    x_solver[6] = pJAK2;
    x_solver[7] = STAT5;
    x_solver[8] = pSTAT5;
    x_solver[9] = SOCS3mRNA;
    x_solver[10] = DecoyR;
    x_solver[11] = IL13_DecoyR;
    x_solver[12] = SOCS3;
    x_solver[13] = CD274mRNA;
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
