# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Isensee_JCB2018.dir/Isensee_JCB2018PYTHON_wrap.cxx"
  "swig/Isensee_JCB2018.py"
  )
endif()
