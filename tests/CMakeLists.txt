cmake_minimum_required (VERSION 2.8)

add_compile_options("-Wall") # caution: might be compiler-specific

#set(CMAKE_CXX_COMPILER g++)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MODULE_PATH ..)

find_package(GSL REQUIRED)
message(STATUS "GSL include dir:" ${GSL_INCLUDE_DIR})
message(STATUS "GSL libraries:" ${GSL_LIBRARIES})
message(STATUS "GSL link directories: " ${GSL_LINK_DIRECTORIES})

include_directories(
    ${GSL_INCLUDE_DIR}
    . .. ../src
)

list(APPEND CMAKE_INSTALL_RPATH ${GSL_LINK_DIRECTORIES})
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

macro(bhm_add_test name src)
  add_executable(${name} ${src})
  add_test(NAME ${name} COMMAND ${name})
  target_link_libraries(${name} samplinglibrary ${GSL_LIBRARIES} libiniparser)
  # install(TARGETS ${name} RUNTIME DESTINATION bin)
endmacro()

# we need a sample input in the test binary directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sample.param ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

bhm_add_test(basic_unittest basic_unittest.cpp)
bhm_add_test(slot_unittest slot_unittest.cpp)
bhm_add_test(histogram_unittest histogram_unittest.cpp)
bhm_add_test(spline_unittest spline_unittest.cpp)
bhm_add_test(iniparser iniparser.cpp)
bhm_add_test(utilities utilities.cpp)
bhm_add_test(read_histogram read_histogram.cpp)
bhm_add_test(logger logger.cpp)
bhm_add_test(grid_test grid_test.cpp)
