# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Oliveira_NatCommun2021.dir/Oliveira_NatCommun2021PYTHON_wrap.cxx"
  "swig/Oliveira_NatCommun2021.py"
  )
endif()
