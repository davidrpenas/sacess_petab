# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "swig/CMakeFiles/_Okuonghae_ChaosSolitonsFractals2020.dir/Okuonghae_ChaosSolitonsFractals2020PYTHON_wrap.cxx"
  "swig/Okuonghae_ChaosSolitonsFractals2020.py"
  )
endif()
