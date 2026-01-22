# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Bertozzi_PNAS2020.py"
  "swig/CMakeFiles/_Bertozzi_PNAS2020.dir/Bertozzi_PNAS2020PYTHON_wrap.cxx"
  )
endif()
