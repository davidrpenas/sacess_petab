#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<std::array<sunindextype, 68>, 67> dJydy_colptrs_Fujita_SciSignal2010_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Fujita_SciSignal2010(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Fujita_SciSignal2010_[index]));
}
} // namespace model_Fujita_SciSignal2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<std::array<sunindextype, 1>, 67> dJydy_rowvals_Fujita_SciSignal2010_ = {{
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
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Fujita_SciSignal2010(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Fujita_SciSignal2010_[index]));
}
} // namespace model_Fujita_SciSignal2010
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
namespace model_Fujita_SciSignal2010 {

void dJydy_Fujita_SciSignal2010(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1 + 1.0*pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3 + 1.0*pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3 + 1.0*pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3 + 1.0*pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0, 2);
            break;
        case 15:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0, 2);
            break;
        case 16:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0, 2);
            break;
        case 18:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0, 2);
            break;
        case 19:
            dJydy[0] = (-1.0*mpAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0 + 1.0*pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0, 2);
            break;
        case 20:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1, 2);
            break;
        case 21:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3, 2);
            break;
        case 22:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0, 2);
            break;
        case 23:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0, 2);
            break;
        case 24:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0, 2);
            break;
        case 25:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0, 2);
            break;
        case 26:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0, 2);
            break;
        case 27:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0, 2);
            break;
        case 28:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0, 2);
            break;
        case 29:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0, 2);
            break;
        case 30:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0, 2);
            break;
        case 31:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0, 2);
            break;
        case 32:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0, 2);
            break;
        case 33:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0, 2);
            break;
        case 34:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0, 2);
            break;
        case 35:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0, 2);
            break;
        case 36:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0, 2);
            break;
        case 37:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0, 2);
            break;
        case 38:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0, 2);
            break;
        case 39:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0, 2);
            break;
        case 40:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0, 2);
            break;
        case 41:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0, 2);
            break;
        case 42:
            dJydy[0] = (-1.0*mpEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0 + 1.0*pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0, 2);
            break;
        case 43:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_01__condition_step_00_1 + 1.0*pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1, 2);
            break;
        case 44:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_01__condition_step_00_3 + 1.0*pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3, 2);
            break;
        case 45:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_01__condition_step_01_0 + 1.0*pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0, 2);
            break;
        case 46:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_01__condition_step_03_0 + 1.0*pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0, 2);
            break;
        case 47:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_01__condition_step_10_0 + 1.0*pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0, 2);
            break;
        case 48:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_01__condition_step_30_0 + 1.0*pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0, 2);
            break;
        case 49:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_02__condition_step_00_1 + 1.0*pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1, 2);
            break;
        case 50:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_02__condition_step_00_3 + 1.0*pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3, 2);
            break;
        case 51:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_02__condition_step_10_0 + 1.0*pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0, 2);
            break;
        case 52:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_02__condition_step_30_0 + 1.0*pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0, 2);
            break;
        case 53:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_025__condition_step_00_3 + 1.0*pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3, 2);
            break;
        case 54:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_025__condition_step_01_0 + 1.0*pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0, 2);
            break;
        case 55:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_03__condition_step_00_3 + 1.0*pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3, 2);
            break;
        case 56:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_03__condition_step_30_0 + 1.0*pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0, 2);
            break;
        case 57:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_035__condition_step_00_1 + 1.0*pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1, 2);
            break;
        case 58:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_035__condition_step_03_0 + 1.0*pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0, 2);
            break;
        case 59:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_035__condition_step_30_0 + 1.0*pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0, 2);
            break;
        case 60:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_04__condition_step_10_0 + 1.0*pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0, 2);
            break;
        case 61:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_07__condition_step_01_0 + 1.0*pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0, 2);
            break;
        case 62:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_07__condition_step_30_0 + 1.0*pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0, 2);
            break;
        case 63:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_075__condition_step_01_0 + 1.0*pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0, 2);
            break;
        case 64:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_08__condition_step_03_0 + 1.0*pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0, 2);
            break;
        case 65:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_08__condition_step_10_0 + 1.0*pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0, 2);
            break;
        case 66:
            dJydy[0] = (-1.0*mpS6_tot__scaling_pS6_tot__0_1__condition_step_00_3 + 1.0*pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3, 2);
            break;
    }
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
