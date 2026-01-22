# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Raimundez_PCB2020.dir/Raimundez_PCB2020PYTHON_wrap.cxx"
  "swig/Raimundez_PCB2020.py"
  )
endif()
