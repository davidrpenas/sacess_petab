# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Zheng_PNAS2012.dir/Zheng_PNAS2012PYTHON_wrap.cxx"
  "swig/Zheng_PNAS2012.py"
  )
endif()
