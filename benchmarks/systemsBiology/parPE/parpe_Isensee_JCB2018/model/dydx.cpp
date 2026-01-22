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
namespace model_Isensee_JCB2018 {

void dydx_Isensee_JCB2018(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[117] = 2*observableParameter1_Calpha_Microscopy*s_Calpha_global;
    dydx[140] = 2*observableParameter1_Calpha_Microscopy*s_Calpha_global;
    dydx[184] = 2*observableParameter1_pRII_Microscopy*rel_open*s_pRII_global;
    dydx[185] = 2*s_pRII_Western;
    dydx[186] = 2*observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global;
    dydx[207] = observableParameter1_pRII_Microscopy*s_pRII_global*(2*rel_open - 2*xi_rel_open*(rel_open - 1));
    dydx[208] = 2*s_pRII_Western;
    dydx[209] = observableParameter1_Calpha_Microscopy*s_Calpha_global*(2*rel_open - 2*xi_rel_open*(rel_open - 1));
    dydx[230] = 2*observableParameter1_pRII_Microscopy*s_pRII_global;
    dydx[231] = 2*s_pRII_Western;
    dydx[345] = 2*observableParameter1_pRII_Microscopy*rel_open*s_pRII_global;
    dydx[346] = 2*s_pRII_Western;
    dydx[347] = 2*observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global;
    dydx[368] = 2*observableParameter1_pRII_Microscopy*rel_open*s_pRII_global;
    dydx[369] = 2*s_pRII_Western;
    dydx[370] = 2*observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global;
    dydx[391] = 2*observableParameter1_pRII_Microscopy*rel_open*s_pRII_global;
    dydx[392] = 2*s_pRII_Western;
    dydx[393] = 2*observableParameter1_Calpha_Microscopy*rel_open*s_Calpha_global;
    dydx[414] = observableParameter1_pRII_Microscopy*s_pRII_global*(2*rel_open - 2*xi_rel_open*(rel_open - 1));
    dydx[415] = 2*s_pRII_Western;
    dydx[416] = observableParameter1_Calpha_Microscopy*s_Calpha_global*(2*rel_open - 2*xi_rel_open*(rel_open - 1));
    dydx[437] = 2*observableParameter1_pRII_Microscopy*s_pRII_global;
    dydx[438] = 2*s_pRII_Western;
    dydx[460] = 2*observableParameter1_pRII_Microscopy*s_pRII_global;
    dydx[461] = 2*s_pRII_Western;
}

} // namespace model_Isensee_JCB2018
} // namespace amici
