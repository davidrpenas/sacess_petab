# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Blasi_CellSystems2016.py"
  "swig/CMakeFiles/_Blasi_CellSystems2016.dir/Blasi_CellSystems2016PYTHON_wrap.cxx"
  )
endif()
