# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Giordano_Nature2020.dir/Giordano_Nature2020PYTHON_wrap.cxx"
  "swig/Giordano_Nature2020.py"
  )
endif()
