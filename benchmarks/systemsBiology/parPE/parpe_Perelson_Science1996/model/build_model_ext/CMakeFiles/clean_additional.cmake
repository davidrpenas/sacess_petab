# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Perelson_Science1996.dir/Perelson_Science1996PYTHON_wrap.cxx"
  "swig/Perelson_Science1996.py"
  )
endif()
