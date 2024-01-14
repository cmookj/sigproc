cmake_minimum_required (VERSION 3.16..3.24)

project (sigproc VERSION 0.0.1
  DESCRIPTION "Signal processing toolkit"
  LANGUAGES CXX)

add_library (libsigproc STATIC
  ${PROJECT_SOURCE_DIR}/src/hampel.cpp)

target_include_directories (libsigproc PUBLIC
  ${PROJECT_SOURCE_DIR}/include)
