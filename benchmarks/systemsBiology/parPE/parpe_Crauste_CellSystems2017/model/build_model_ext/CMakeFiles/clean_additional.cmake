# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Crauste_CellSystems2017.dir/Crauste_CellSystems2017PYTHON_wrap.cxx"
  "swig/Crauste_CellSystems2017.py"
  )
endif()
