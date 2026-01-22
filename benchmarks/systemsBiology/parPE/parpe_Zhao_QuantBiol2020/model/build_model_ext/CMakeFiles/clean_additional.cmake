# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Zhao_QuantBiol2020.dir/Zhao_QuantBiol2020PYTHON_wrap.cxx"
  "swig/Zhao_QuantBiol2020.py"
  )
endif()
