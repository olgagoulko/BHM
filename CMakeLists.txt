cmake_minimum_required (VERSION 2.8)
project(sampling)

option(ENABLE_TESTING "Build and enable unit tests" ON)
mark_as_advanced(ENABLE_TESTING) # users normally should always enable testing!

add_subdirectory(iniparser)
add_subdirectory (src)
if (ENABLE_TESTING)
  enable_testing()
  add_subdirectory (tests)
endif()

#add_test(test basic_test)
