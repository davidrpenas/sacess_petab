# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Rahman_MBS2016.dir/Rahman_MBS2016PYTHON_wrap.cxx"
  "swig/Rahman_MBS2016.py"
  )
endif()
