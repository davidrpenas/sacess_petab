# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/Armistead_CellDeathDis2024.py"
  "swig/CMakeFiles/_Armistead_CellDeathDis2024.dir/Armistead_CellDeathDis2024PYTHON_wrap.cxx"
  )
endif()
