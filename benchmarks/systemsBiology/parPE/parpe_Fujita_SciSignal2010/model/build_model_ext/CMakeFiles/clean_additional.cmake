# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Fujita_SciSignal2010.dir/Fujita_SciSignal2010PYTHON_wrap.cxx"
  "swig/Fujita_SciSignal2010.py"
  )
endif()
