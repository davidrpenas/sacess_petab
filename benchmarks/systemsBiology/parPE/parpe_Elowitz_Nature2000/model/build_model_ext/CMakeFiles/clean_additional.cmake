# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Elowitz_Nature2000.dir/Elowitz_Nature2000PYTHON_wrap.cxx"
  "swig/Elowitz_Nature2000.py"
  )
endif()
