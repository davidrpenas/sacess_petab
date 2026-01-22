# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Boehm_JProteomeRes2014.py"
  "swig/CMakeFiles/_Boehm_JProteomeRes2014.dir/Boehm_JProteomeRes2014PYTHON_wrap.cxx"
  )
endif()
