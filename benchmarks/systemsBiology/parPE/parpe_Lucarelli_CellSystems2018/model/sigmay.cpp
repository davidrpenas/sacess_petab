#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void sigmay_Lucarelli_CellSystems2018(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_Bmp4_gene_expression_OE_nExpID100 = sd_Bmp4_nExpID100;  // sigmay[0]
    sigma_Cxcl15_gene_expression_OE_nExpID100 = sd_Cxcl15_nExpID100;  // sigmay[1]
    sigma_Dnmt3a_gene_expression_OE_nExpID100 = sd_Dnmt3a_nExpID100;  // sigmay[2]
    sigma_Dusp5_gene_expression_OE_nExpID100 = sd_Dusp5_nExpID100;  // sigmay[3]
    sigma_Jun_gene_expression_OE_nExpID100 = sd_Jun_nExpID100;  // sigmay[4]
    sigma_Klf10_gene_expression_OE_nExpID100 = sd_Klf10_nExpID100;  // sigmay[5]
    sigma_Pdk4_gene_expression_OE_nExpID100 = sd_Pdk4_nExpID100;  // sigmay[6]
    sigma_S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 = 5.0;  // sigmay[7]
    sigma_S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[8]
    sigma_S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 = 5.0;  // sigmay[9]
    sigma_S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[10]
    sigma_S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 = 5.0;  // sigmay[11]
    sigma_S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[12]
    sigma_S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 = 5.0;  // sigmay[13]
    sigma_S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[14]
    sigma_S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 = 5.0;  // sigmay[15]
    sigma_S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[16]
    sigma_S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 = 5.0;  // sigmay[17]
    sigma_S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[18]
    sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 = 15.0;  // sigmay[19]
    sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 = 1.0;  // sigmay[20]
    sigma_S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[21]
    sigma_S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 = 1.0;  // sigmay[22]
    sigma_S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[23]
    sigma_S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[24]
    sigma_S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[25]
    sigma_S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[26]
    sigma_S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[27]
    sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[28]
    sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[29]
    sigma_S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[30]
    sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[31]
    sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[32]
    sigma_S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[33]
    sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[34]
    sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[35]
    sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[36]
    sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[37]
    sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[38]
    sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[39]
    sigma_S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[40]
    sigma_S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[41]
    sigma_S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[42]
    sigma_S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[43]
    sigma_S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[44]
    sigma_S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[45]
    sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[46]
    sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[47]
    sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[48]
    sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[49]
    sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[50]
    sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[51]
    sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[52]
    sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[53]
    sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[54]
    sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[55]
    sigma_S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 = 15.0;  // sigmay[56]
    sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 = 5.0;  // sigmay[57]
    sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 = 5.0;  // sigmay[58]
    sigma_Ski_gene_expression_OE_nExpID100 = sd_Ski_nExpID100;  // sigmay[59]
    sigma_Skil_gene_expression_OE_nExpID100 = sd_Skil_nExpID100;  // sigmay[60]
    sigma_Smad7_gene_expression_OE_nExpID100 = sd_Smad7_nExpID100;  // sigmay[61]
    sigma_Sox4_gene_expression_OE_nExpID100 = sd_Sox4_nExpID100;  // sigmay[62]
    sigma_TGFbR_PL_PN_TGFbR_hepa16_nExpID11 = 0.091999999999999998;  // sigmay[63]
    sigma_Tgfa_gene_expression_OE_nExpID100 = sd_Tgfa_nExpID100;  // sigmay[64]
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
