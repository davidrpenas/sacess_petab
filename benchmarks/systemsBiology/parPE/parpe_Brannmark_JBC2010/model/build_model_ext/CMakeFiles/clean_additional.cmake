# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Brannmark_JBC2010.py"
  "swig/CMakeFiles/_Brannmark_JBC2010.dir/Brannmark_JBC2010PYTHON_wrap.cxx"
  )
endif()
