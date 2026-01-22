#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void dJydsigma_Lucarelli_CellSystems2018(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_Bmp4_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Bmp4_gene_expression_OE_nExpID100)/M_LN10 - std::log(mBmp4_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Bmp4_gene_expression_OE_nExpID100, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_Cxcl15_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Cxcl15_gene_expression_OE_nExpID100)/M_LN10 - std::log(mCxcl15_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Cxcl15_gene_expression_OE_nExpID100, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_Dnmt3a_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Dnmt3a_gene_expression_OE_nExpID100)/M_LN10 - std::log(mDnmt3a_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Dnmt3a_gene_expression_OE_nExpID100, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_Dusp5_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Dusp5_gene_expression_OE_nExpID100)/M_LN10 - std::log(mDusp5_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Dusp5_gene_expression_OE_nExpID100, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_Jun_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Jun_gene_expression_OE_nExpID100)/M_LN10 - std::log(mJun_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Jun_gene_expression_OE_nExpID100, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_Klf10_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Klf10_gene_expression_OE_nExpID100)/M_LN10 - std::log(mKlf10_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Klf10_gene_expression_OE_nExpID100, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_Pdk4_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Pdk4_gene_expression_OE_nExpID100)/M_LN10 - std::log(mPdk4_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Pdk4_gene_expression_OE_nExpID100, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*std::pow(S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*std::pow(S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*std::pow(S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*std::pow(S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*std::pow(S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*std::pow(S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - 1.0*std::pow(S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - mS23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2)/std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - 1.0*std::pow(S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - mS23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2)/std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - 1.0*std::pow(S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - mS23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2)/std::pow(sigma_S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - 1.0*std::pow(S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - mS4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2)/std::pow(sigma_S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*std::pow(S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*std::pow(S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_Ski_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Ski_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSki_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Ski_gene_expression_OE_nExpID100, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_Skil_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Skil_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSkil_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Skil_gene_expression_OE_nExpID100, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_Smad7_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Smad7_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSmad7_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Smad7_gene_expression_OE_nExpID100, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_Sox4_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Sox4_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSox4_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Sox4_gene_expression_OE_nExpID100, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_TGFbR_PL_PN_TGFbR_hepa16_nExpID11 - 1.0*std::pow(TGFbR_PL_PN_TGFbR_hepa16_nExpID11 - mTGFbR_PL_PN_TGFbR_hepa16_nExpID11, 2)/std::pow(sigma_TGFbR_PL_PN_TGFbR_hepa16_nExpID11, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_Tgfa_gene_expression_OE_nExpID100 - 1.0*std::pow(std::log(Tgfa_gene_expression_OE_nExpID100)/M_LN10 - std::log(mTgfa_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Tgfa_gene_expression_OE_nExpID100, 3);
            break;
    }
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
