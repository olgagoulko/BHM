cmake_minimum_required (VERSION 2.8)

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MODULE_PATH ..)

find_package(GSL REQUIRED)
message("GSL include dir:" ${GSL_INCLUDE_DIR})
message("GSL libraries:" ${GSL_LIBRARIES})

include_directories(
    ${GSL_INCLUDE_DIR}
    . ..
    )
    
add_library (samplinglibrary basic.cpp matrix.cpp slot.cpp histogram.cpp histogram.hpp slot.hpp matrix.hpp basic.hpp)
add_executable(main main.cpp)
    
target_link_libraries(main samplinglibrary ${GSL_LIBRARIES})

install(TARGETS main RUNTIME DESTINATION bin)
