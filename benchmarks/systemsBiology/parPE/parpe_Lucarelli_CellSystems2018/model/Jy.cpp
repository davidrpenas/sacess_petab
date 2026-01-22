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

void Jy_Lucarelli_CellSystems2018(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mBmp4_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Bmp4_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Bmp4_gene_expression_OE_nExpID100)/M_LN10 - std::log(mBmp4_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Bmp4_gene_expression_OE_nExpID100, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mCxcl15_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Cxcl15_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Cxcl15_gene_expression_OE_nExpID100)/M_LN10 - std::log(mCxcl15_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Cxcl15_gene_expression_OE_nExpID100, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mDnmt3a_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Dnmt3a_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Dnmt3a_gene_expression_OE_nExpID100)/M_LN10 - std::log(mDnmt3a_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Dnmt3a_gene_expression_OE_nExpID100, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mDusp5_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Dusp5_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Dusp5_gene_expression_OE_nExpID100)/M_LN10 - std::log(mDusp5_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Dusp5_gene_expression_OE_nExpID100, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mJun_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Jun_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Jun_gene_expression_OE_nExpID100)/M_LN10 - std::log(mJun_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Jun_gene_expression_OE_nExpID100, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mKlf10_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Klf10_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Klf10_gene_expression_OE_nExpID100)/M_LN10 - std::log(mKlf10_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Klf10_gene_expression_OE_nExpID100, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mPdk4_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Pdk4_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Pdk4_gene_expression_OE_nExpID100)/M_LN10 - std::log(mPdk4_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Pdk4_gene_expression_OE_nExpID100, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)) + 0.5*std::pow(S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)) + 0.5*std::pow(S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)) + 0.5*std::pow(S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)) + 0.5*std::pow(S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)) + 0.5*std::pow(S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)) + 0.5*std::pow(S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - mS23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2)/std::pow(sigma_S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2)) + 0.5*std::pow(S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - mS23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2)/std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2)) + 0.5*std::pow(S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - mS23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2)/std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2)) + 0.5*std::pow(S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - mS23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2)/std::pow(sigma_S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2)) + 0.5*std::pow(S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - mS4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2)/std::pow(sigma_S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)) + 0.5*std::pow(S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - mS4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2)/std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)) + 0.5*std::pow(S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - mS4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2)/std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mSki_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Ski_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Ski_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSki_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Ski_gene_expression_OE_nExpID100, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mSkil_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Skil_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Skil_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSkil_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Skil_gene_expression_OE_nExpID100, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mSmad7_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Smad7_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Smad7_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSmad7_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Smad7_gene_expression_OE_nExpID100, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mSox4_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Sox4_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Sox4_gene_expression_OE_nExpID100)/M_LN10 - std::log(mSox4_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Sox4_gene_expression_OE_nExpID100, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_TGFbR_PL_PN_TGFbR_hepa16_nExpID11, 2)) + 0.5*std::pow(TGFbR_PL_PN_TGFbR_hepa16_nExpID11 - mTGFbR_PL_PN_TGFbR_hepa16_nExpID11, 2)/std::pow(sigma_TGFbR_PL_PN_TGFbR_hepa16_nExpID11, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mTgfa_gene_expression_OE_nExpID100, 2)*std::pow(sigma_Tgfa_gene_expression_OE_nExpID100, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Tgfa_gene_expression_OE_nExpID100)/M_LN10 - std::log(mTgfa_gene_expression_OE_nExpID100)/M_LN10, 2)/std::pow(sigma_Tgfa_gene_expression_OE_nExpID100, 2);
            break;
    }
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
