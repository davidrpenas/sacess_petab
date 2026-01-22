# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Fiedler_BMCSystBiol2016.dir/Fiedler_BMCSystBiol2016PYTHON_wrap.cxx"
  "swig/Fiedler_BMCSystBiol2016.py"
  )
endif()
