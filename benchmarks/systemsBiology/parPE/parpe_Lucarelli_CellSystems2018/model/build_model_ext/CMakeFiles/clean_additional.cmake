# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Lucarelli_CellSystems2018.dir/Lucarelli_CellSystems2018PYTHON_wrap.cxx"
  "swig/Lucarelli_CellSystems2018.py"
  )
endif()
