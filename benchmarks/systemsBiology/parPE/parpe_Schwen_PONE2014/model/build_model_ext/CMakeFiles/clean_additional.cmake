# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Schwen_PONE2014.dir/Schwen_PONE2014PYTHON_wrap.cxx"
  "swig/Schwen_PONE2014.py"
  )
endif()
