# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Smith_BMCSystBiol2013.dir/Smith_BMCSystBiol2013PYTHON_wrap.cxx"
  "swig/Smith_BMCSystBiol2013.py"
  )
endif()
