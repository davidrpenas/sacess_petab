# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Raia_CancerResearch2011.dir/Raia_CancerResearch2011PYTHON_wrap.cxx"
  "swig/Raia_CancerResearch2011.py"
  )
endif()
