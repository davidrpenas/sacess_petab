# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Bruno_JExpBot2016.py"
  "swig/CMakeFiles/_Bruno_JExpBot2016.dir/Bruno_JExpBot2016PYTHON_wrap.cxx"
  )
endif()
