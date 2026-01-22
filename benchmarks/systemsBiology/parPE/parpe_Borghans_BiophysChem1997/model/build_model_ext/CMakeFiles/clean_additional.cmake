# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Borghans_BiophysChem1997.py"
  "swig/CMakeFiles/_Borghans_BiophysChem1997.dir/Borghans_BiophysChem1997PYTHON_wrap.cxx"
  )
endif()
