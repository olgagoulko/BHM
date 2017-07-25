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
message("GSL include dir:" ${GSL_INCLUDE_DIR})
message("GSL libraries:" ${GSL_LIBRARIES})

include_directories(
    ${GSL_INCLUDE_DIR}
    . .. ../src
    )

add_executable(maintest basic_unittest.cpp sput.h)
add_test(NAME maintest COMMAND maintest)

target_link_libraries(maintest samplinglibrary ${GSL_LIBRARIES})

install(TARGETS maintest RUNTIME DESTINATION bin)
