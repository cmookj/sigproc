cmake_minimum_required (VERSION 3.16..3.24)

project (signal_processing VERSION 0.0.1
  DESCRIPTION "Signal processing toolkit"
  LANGUAGES CXX)

add_library (libsignal STATIC
  ${PROJECT_SOURCE_DIR}/src/hampel.cpp)

target_include_directories (libsignal PUBLIC
  ${PROJECT_SOURCE_DIR}/include)
