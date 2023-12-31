cmake_minimum_required (VERSION 3.16..3.24)

project (hampel VERSION 0.0.1
  DESCRIPTION "Hampel filter"
  LANGUAGES CXX)

add_library (libhampel STATIC
  ${PROJECT_SOURCE_DIR}/src/hampel.cpp)

target_include_directories (libhampel PUBLIC 
  ${PROJECT_SOURCE_DIR}/include)
