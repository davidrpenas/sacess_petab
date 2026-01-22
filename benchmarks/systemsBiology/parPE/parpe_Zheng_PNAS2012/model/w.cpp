#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void w_Zheng_PNAS2012(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        inflow = dilution*inflowp;  // w[0]
        flux_v1_v_0 = inflow;  // w[1]
    }

    // dynamic expressions
    flux_v2_v_1 = K27me0K36me0*k00_01;  // w[2]
    flux_v3_v_2 = K27me0K36me1*k01_00;  // w[3]
    flux_v4_v_3 = K27me0K36me0*k00_10;  // w[4]
    flux_v5_v_4 = K27me1K36me0*k10_00;  // w[5]
    flux_v6_v_5 = K27me1K36me0*k10_11;  // w[6]
    flux_v7_v_6 = K27me1K36me1*k11_10;  // w[7]
    flux_v8_v_7 = K27me0K36me1*k01_11;  // w[8]
    flux_v9_v_8 = K27me1K36me1*k11_01;  // w[9]
    flux_v10_v_9 = K27me0K36me1*k01_02;  // w[10]
    flux_v11_v_10 = K27me0K36me2*k02_01;  // w[11]
    flux_v12_v_11 = K27me1K36me0*k10_20;  // w[12]
    flux_v13_v_12 = K27me2K36me0*k20_10;  // w[13]
    flux_v14_v_13 = K27me0K36me2*k02_12;  // w[14]
    flux_v15_v_14 = K27me1K36me2*k12_02;  // w[15]
    flux_v16_v_15 = K27me2K36me0*k20_21;  // w[16]
    flux_v17_v_16 = K27me2K36me1*k21_20;  // w[17]
    flux_v18_v_17 = K27me0K36me2*k02_03;  // w[18]
    flux_v19_v_18 = K27me0K36me3*k03_02;  // w[19]
    flux_v20_v_19 = K27me2K36me0*k20_30;  // w[20]
    flux_v21_v_20 = K27me3K36me0*k30_20;  // w[21]
    flux_v22_v_21 = K27me1K36me1*k11_12;  // w[22]
    flux_v23_v_22 = K27me1K36me2*k12_11;  // w[23]
    flux_v24_v_23 = K27me1K36me1*k11_21;  // w[24]
    flux_v25_v_24 = K27me2K36me1*k21_11;  // w[25]
    flux_v26_v_25 = K27me1K36me2*k12_22;  // w[26]
    flux_v27_v_26 = K27me2K36me2*k22_12;  // w[27]
    flux_v28_v_27 = K27me2K36me1*k21_22;  // w[28]
    flux_v29_v_28 = K27me2K36me2*k22_21;  // w[29]
    flux_v30_v_29 = K27me2K36me2*k22_23;  // w[30]
    flux_v31_v_30 = K27me2K36me3*k23_22;  // w[31]
    flux_v32_v_31 = K27me2K36me2*k22_32;  // w[32]
    flux_v33_v_32 = K27me3K36me2*k32_22;  // w[33]
    flux_v34_v_33 = K27me3K36me0*k30_31;  // w[34]
    flux_v35_v_34 = K27me3K36me1*k31_30;  // w[35]
    flux_v36_v_35 = K27me0K36me3*k03_13;  // w[36]
    flux_v37_v_36 = K27me1K36me3*k13_03;  // w[37]
    flux_v38_v_37 = K27me1K36me2*k12_13;  // w[38]
    flux_v39_v_38 = K27me1K36me3*k13_12;  // w[39]
    flux_v40_v_39 = K27me2K36me1*k21_31;  // w[40]
    flux_v41_v_40 = K27me3K36me1*k31_21;  // w[41]
    flux_v42_v_41 = K27me3K36me1*k31_32;  // w[42]
    flux_v43_v_42 = K27me3K36me2*k32_31;  // w[43]
    flux_v44_v_43 = K27me1K36me3*k13_23;  // w[44]
    flux_v45_v_44 = K27me2K36me3*k23_13;  // w[45]
    flux_v46_v_45 = K27me0K36me0*inflow;  // w[46]
    flux_v47_v_46 = K27me0K36me1*inflow;  // w[47]
    flux_v48_v_47 = K27me1K36me0*inflow;  // w[48]
    flux_v49_v_48 = K27me0K36me2*inflow;  // w[49]
    flux_v50_v_49 = K27me1K36me1*inflow;  // w[50]
    flux_v51_v_50 = K27me2K36me0*inflow;  // w[51]
    flux_v52_v_51 = K27me0K36me3*inflow;  // w[52]
    flux_v53_v_52 = K27me1K36me2*inflow;  // w[53]
    flux_v54_v_53 = K27me2K36me1*inflow;  // w[54]
    flux_v55_v_54 = K27me3K36me0*inflow;  // w[55]
    flux_v56_v_55 = K27me1K36me3*inflow;  // w[56]
    flux_v57_v_56 = K27me2K36me2*inflow;  // w[57]
    flux_v58_v_57 = K27me3K36me1*inflow;  // w[58]
    flux_v59_v_58 = K27me2K36me3*inflow;  // w[59]
    flux_v60_v_59 = K27me3K36me2*inflow;  // w[60]
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
