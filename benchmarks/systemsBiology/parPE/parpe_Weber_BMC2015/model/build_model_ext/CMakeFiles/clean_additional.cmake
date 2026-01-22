# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Weber_BMC2015.dir/Weber_BMC2015PYTHON_wrap.cxx"
  "swig/Weber_BMC2015.py"
  )
endif()
