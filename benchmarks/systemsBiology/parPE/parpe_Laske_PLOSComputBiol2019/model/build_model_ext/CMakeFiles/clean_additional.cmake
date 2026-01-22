# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Laske_PLOSComputBiol2019.dir/Laske_PLOSComputBiol2019PYTHON_wrap.cxx"
  "swig/Laske_PLOSComputBiol2019.py"
  )
endif()
