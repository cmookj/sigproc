# BUILD_TESTING variable is created by include (CTest)
# It is set to ON by default 

if (BUILD_TESTING)
  add_executable (hampel_test hampel_test.cpp)
  target_link_libraries (hampel_test PRIVATE 
    libhampel 
    /usr/local/lib/librbm.a
    /usr/local/lib/libblat.a
    GTest::gtest_main)

  target_include_directories (hampel_test PRIVATE 
    ${PROJECT_SOURCE_DIR}/hampel/include
    /usr/local/include
    )

  include (GoogleTest)

  # Finds all the Google tests associated with the executable 
  gtest_discover_tests (hampel_test)

endif ()
