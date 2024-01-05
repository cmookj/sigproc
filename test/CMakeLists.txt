# BUILD_TESTING variable is created by include (CTest)
# It is set to ON by default

if (BUILD_TESTING)
  add_executable (signal_processing_test signal_processing_test.cpp)
  target_link_libraries (signal_processing_test PRIVATE
    libsignal
    /usr/local/lib/librbm.a
    /usr/local/lib/libblat.a
    GTest::gtest_main)

  target_include_directories (signal_processing_test PRIVATE
    ${PROJECT_SOURCE_DIR}/signal_processing/include
    /usr/local/include
    )

  include (GoogleTest)

  # Finds all the Google tests associated with the executable
  gtest_discover_tests (signal_processing_test)

endif ()
