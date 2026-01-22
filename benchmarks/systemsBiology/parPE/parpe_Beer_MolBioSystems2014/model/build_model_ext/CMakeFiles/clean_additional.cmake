# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Beer_MolBioSystems2014.py"
  "swig/CMakeFiles/_Beer_MolBioSystems2014.dir/Beer_MolBioSystems2014PYTHON_wrap.cxx"
  )
endif()
