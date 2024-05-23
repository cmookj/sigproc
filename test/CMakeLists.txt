# BUILD_TESTING variable is created by include (CTest)
# It is set to ON by default

if (BUILD_TESTING)
  add_executable (sigproc_test
    hampel_test.cpp
    tail_window_stat_test.cpp
  )
  target_link_libraries (sigproc_test PRIVATE
    libsigproc
    GTest::gtest_main)

  target_include_directories (sigproc_test PRIVATE
    ${PROJECT_SOURCE_DIR}/sigproc/include
    /usr/local/include
    )

  include (GoogleTest)

  # Finds all the Google tests associated with the executable
  gtest_discover_tests (sigproc_test)

endif ()
