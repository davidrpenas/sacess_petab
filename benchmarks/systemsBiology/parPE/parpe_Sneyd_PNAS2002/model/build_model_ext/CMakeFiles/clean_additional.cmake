# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Sneyd_PNAS2002.dir/Sneyd_PNAS2002PYTHON_wrap.cxx"
  "swig/Sneyd_PNAS2002.py"
  )
endif()
