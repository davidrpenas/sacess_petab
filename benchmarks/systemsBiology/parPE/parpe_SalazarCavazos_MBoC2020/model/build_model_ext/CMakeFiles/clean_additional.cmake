# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_SalazarCavazos_MBoC2020.dir/SalazarCavazos_MBoC2020PYTHON_wrap.cxx"
  "swig/SalazarCavazos_MBoC2020.py"
  )
endif()
