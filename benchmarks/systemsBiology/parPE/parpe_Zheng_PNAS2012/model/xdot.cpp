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
namespace model_Zheng_PNAS2012 {

void xdot_Zheng_PNAS2012(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dK27me0K36me0dt = flux_v1_v_0 - flux_v2_v_1 + flux_v3_v_2 - flux_v46_v_45 - flux_v4_v_3 + flux_v5_v_4;  // xdot[0]
    dK27me0K36me1dt = -flux_v10_v_9 + flux_v11_v_10 + flux_v2_v_1 - flux_v3_v_2 - flux_v47_v_46 - flux_v8_v_7 + flux_v9_v_8;  // xdot[1]
    dK27me1K36me0dt = -flux_v12_v_11 + flux_v13_v_12 - flux_v48_v_47 + flux_v4_v_3 - flux_v5_v_4 - flux_v6_v_5 + flux_v7_v_6;  // xdot[2]
    dK27me0K36me2dt = flux_v10_v_9 - flux_v11_v_10 - flux_v14_v_13 + flux_v15_v_14 - flux_v18_v_17 + flux_v19_v_18 - flux_v49_v_48;  // xdot[3]
    dK27me1K36me1dt = -flux_v22_v_21 + flux_v23_v_22 - flux_v24_v_23 + flux_v25_v_24 - flux_v50_v_49 + flux_v6_v_5 - flux_v7_v_6 + flux_v8_v_7 - flux_v9_v_8;  // xdot[4]
    dK27me2K36me0dt = flux_v12_v_11 - flux_v13_v_12 - flux_v16_v_15 + flux_v17_v_16 - flux_v20_v_19 + flux_v21_v_20 - flux_v51_v_50;  // xdot[5]
    dK27me0K36me3dt = flux_v18_v_17 - flux_v19_v_18 - flux_v36_v_35 + flux_v37_v_36 - flux_v52_v_51;  // xdot[6]
    dK27me1K36me2dt = flux_v14_v_13 - flux_v15_v_14 + flux_v22_v_21 - flux_v23_v_22 - flux_v26_v_25 + flux_v27_v_26 - flux_v38_v_37 + flux_v39_v_38 - flux_v53_v_52;  // xdot[7]
    dK27me2K36me1dt = flux_v16_v_15 - flux_v17_v_16 + flux_v24_v_23 - flux_v25_v_24 - flux_v28_v_27 + flux_v29_v_28 - flux_v40_v_39 + flux_v41_v_40 - flux_v54_v_53;  // xdot[8]
    dK27me3K36me0dt = flux_v20_v_19 - flux_v21_v_20 - flux_v34_v_33 + flux_v35_v_34 - flux_v55_v_54;  // xdot[9]
    dK27me1K36me3dt = flux_v36_v_35 - flux_v37_v_36 + flux_v38_v_37 - flux_v39_v_38 - flux_v44_v_43 + flux_v45_v_44 - flux_v56_v_55;  // xdot[10]
    dK27me2K36me2dt = flux_v26_v_25 - flux_v27_v_26 + flux_v28_v_27 - flux_v29_v_28 - flux_v30_v_29 + flux_v31_v_30 - flux_v32_v_31 + flux_v33_v_32 - flux_v57_v_56;  // xdot[11]
    dK27me3K36me1dt = flux_v34_v_33 - flux_v35_v_34 + flux_v40_v_39 - flux_v41_v_40 - flux_v42_v_41 + flux_v43_v_42 - flux_v58_v_57;  // xdot[12]
    dK27me2K36me3dt = flux_v30_v_29 - flux_v31_v_30 + flux_v44_v_43 - flux_v45_v_44 - flux_v59_v_58;  // xdot[13]
    dK27me3K36me2dt = flux_v32_v_31 - flux_v33_v_32 + flux_v42_v_41 - flux_v43_v_42 - flux_v60_v_59;  // xdot[14]
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
