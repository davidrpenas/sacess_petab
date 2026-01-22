#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<std::array<sunindextype, 66>, 65> dJydy_colptrs_Lucarelli_CellSystems2018_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Lucarelli_CellSystems2018(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Lucarelli_CellSystems2018_[index]));
}
} // namespace model_Lucarelli_CellSystems2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<std::array<sunindextype, 1>, 65> dJydy_rowvals_Lucarelli_CellSystems2018_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Lucarelli_CellSystems2018(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Lucarelli_CellSystems2018_[index]));
}
} // namespace model_Lucarelli_CellSystems2018
} // namespace amici




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
#include "dJydy.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void dJydy_Lucarelli_CellSystems2018(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*std::log(Bmp4_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mBmp4_gene_expression_OE_nExpID100)/M_LN10)/(Bmp4_gene_expression_OE_nExpID100*std::pow(sigma_Bmp4_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 1:
            dJydy[0] = (1.0*std::log(Cxcl15_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mCxcl15_gene_expression_OE_nExpID100)/M_LN10)/(Cxcl15_gene_expression_OE_nExpID100*std::pow(sigma_Cxcl15_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 2:
            dJydy[0] = (1.0*std::log(Dnmt3a_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mDnmt3a_gene_expression_OE_nExpID100)/M_LN10)/(Dnmt3a_gene_expression_OE_nExpID100*std::pow(sigma_Dnmt3a_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 3:
            dJydy[0] = (1.0*std::log(Dusp5_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mDusp5_gene_expression_OE_nExpID100)/M_LN10)/(Dusp5_gene_expression_OE_nExpID100*std::pow(sigma_Dusp5_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 4:
            dJydy[0] = (1.0*std::log(Jun_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mJun_gene_expression_OE_nExpID100)/M_LN10)/(Jun_gene_expression_OE_nExpID100*std::pow(sigma_Jun_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 5:
            dJydy[0] = (1.0*std::log(Klf10_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mKlf10_gene_expression_OE_nExpID100)/M_LN10)/(Klf10_gene_expression_OE_nExpID100*std::pow(sigma_Klf10_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 6:
            dJydy[0] = (1.0*std::log(Pdk4_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mPdk4_gene_expression_OE_nExpID100)/M_LN10)/(Pdk4_gene_expression_OE_nExpID100*std::pow(sigma_Pdk4_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 7:
            dJydy[0] = (1.0*S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*mS23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5)/std::pow(sigma_S23IP_S2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 8:
            dJydy[0] = (1.0*S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 9:
            dJydy[0] = (1.0*S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*mS23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5)/std::pow(sigma_S23IP_S3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 10:
            dJydy[0] = (1.0*S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 11:
            dJydy[0] = (1.0*S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*mS23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5)/std::pow(sigma_S23IP_pS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 12:
            dJydy[0] = (1.0*S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 13:
            dJydy[0] = (1.0*S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*mS23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5)/std::pow(sigma_S23IP_pS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 14:
            dJydy[0] = (1.0*S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 15:
            dJydy[0] = (1.0*S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*mS23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5)/std::pow(sigma_S23IP_ppS2_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 16:
            dJydy[0] = (1.0*S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 17:
            dJydy[0] = (1.0*S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5 - 1.0*mS23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5)/std::pow(sigma_S23IP_ppS3_MSPL2_PL_hepa16_pS2_pS3_time_course_nExpID5, 2);
            break;
        case 18:
            dJydy[0] = (1.0*S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 19:
            dJydy[0] = (1.0*S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - 1.0*mS23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1)/std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2);
            break;
        case 20:
            dJydy[0] = (1.0*S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - 1.0*mS23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10)/std::pow(sigma_S23IP_tS2_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2);
            break;
        case 21:
            dJydy[0] = (1.0*S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 22:
            dJydy[0] = (1.0*S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10 - 1.0*mS23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10)/std::pow(sigma_S23IP_tS3_MSPL2_PL_PN_molecules_per_cell_hepa16_ratios_nExpID10, 2);
            break;
        case 23:
            dJydy[0] = (1.0*S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 24:
            dJydy[0] = (1.0*S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S23IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 25:
            dJydy[0] = (1.0*S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 26:
            dJydy[0] = (1.0*S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S2IP_S2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 27:
            dJydy[0] = (1.0*S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 28:
            dJydy[0] = (1.0*S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 29:
            dJydy[0] = (1.0*S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S2IP_pS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 30:
            dJydy[0] = (1.0*S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 31:
            dJydy[0] = (1.0*S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 32:
            dJydy[0] = (1.0*S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S2IP_ppS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 33:
            dJydy[0] = (1.0*S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 34:
            dJydy[0] = (1.0*S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 35:
            dJydy[0] = (1.0*S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S2IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 36:
            dJydy[0] = (1.0*S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 37:
            dJydy[0] = (1.0*S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S2IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 38:
            dJydy[0] = (1.0*S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 39:
            dJydy[0] = (1.0*S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S2IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 40:
            dJydy[0] = (1.0*S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_S2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 41:
            dJydy[0] = (1.0*S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_S3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 42:
            dJydy[0] = (1.0*S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_pS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 43:
            dJydy[0] = (1.0*S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_pS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 44:
            dJydy[0] = (1.0*S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_ppS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 45:
            dJydy[0] = (1.0*S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_ppS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 46:
            dJydy[0] = (1.0*S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 47:
            dJydy[0] = (1.0*S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S3IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 48:
            dJydy[0] = (1.0*S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 49:
            dJydy[0] = (1.0*S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S3IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 50:
            dJydy[0] = (1.0*S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 51:
            dJydy[0] = (1.0*S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S3IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 52:
            dJydy[0] = (1.0*S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 53:
            dJydy[0] = (1.0*S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S4IP_tS2_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 54:
            dJydy[0] = (1.0*S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 55:
            dJydy[0] = (1.0*S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S4IP_tS3_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 56:
            dJydy[0] = (1.0*S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1 - 1.0*mS4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1)/std::pow(sigma_S4IP_tS4_MSPL2_PL_PN_molecules_per_cell_hepa16_nExpID1, 2);
            break;
        case 57:
            dJydy[0] = (1.0*S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2 - 1.0*mS4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2)/std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_no_single_complex_nExpID2, 2);
            break;
        case 58:
            dJydy[0] = (1.0*S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6 - 1.0*mS4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6)/std::pow(sigma_S4IP_tS4_MSPL2_complexes_all_hepa16_time_course_nExpID6, 2);
            break;
        case 59:
            dJydy[0] = (1.0*std::log(Ski_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mSki_gene_expression_OE_nExpID100)/M_LN10)/(Ski_gene_expression_OE_nExpID100*std::pow(sigma_Ski_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 60:
            dJydy[0] = (1.0*std::log(Skil_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mSkil_gene_expression_OE_nExpID100)/M_LN10)/(Skil_gene_expression_OE_nExpID100*std::pow(sigma_Skil_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 61:
            dJydy[0] = (1.0*std::log(Smad7_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mSmad7_gene_expression_OE_nExpID100)/M_LN10)/(Smad7_gene_expression_OE_nExpID100*std::pow(sigma_Smad7_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 62:
            dJydy[0] = (1.0*std::log(Sox4_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mSox4_gene_expression_OE_nExpID100)/M_LN10)/(Sox4_gene_expression_OE_nExpID100*std::pow(sigma_Sox4_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
        case 63:
            dJydy[0] = (1.0*TGFbR_PL_PN_TGFbR_hepa16_nExpID11 - 1.0*mTGFbR_PL_PN_TGFbR_hepa16_nExpID11)/std::pow(sigma_TGFbR_PL_PN_TGFbR_hepa16_nExpID11, 2);
            break;
        case 64:
            dJydy[0] = (1.0*std::log(Tgfa_gene_expression_OE_nExpID100)/M_LN10 - 1.0*std::log(mTgfa_gene_expression_OE_nExpID100)/M_LN10)/(Tgfa_gene_expression_OE_nExpID100*std::pow(sigma_Tgfa_gene_expression_OE_nExpID100, 2)*M_LN10);
            break;
    }
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
