cmake_minimum_required (VERSION 2.6)
project(sampling)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(sampling main.cpp basic.cpp)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(GSL REQUIRED)
message("GSL include dir:" ${GSL_INCLUDE_DIR})
message("GSL libraries:" ${GSL_LIBRARIES})

include_directories(
    ${GSL_INCLUDE_DIR}
    )
target_link_libraries(sampling ${GSL_LIBRARIES})

install(TARGETS sampling RUNTIME DESTINATION bin)
