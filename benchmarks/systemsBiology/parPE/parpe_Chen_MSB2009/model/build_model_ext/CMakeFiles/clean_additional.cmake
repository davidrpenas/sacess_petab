# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Chen_MSB2009.dir/Chen_MSB2009PYTHON_wrap.cxx"
  "swig/Chen_MSB2009.py"
  )
endif()
