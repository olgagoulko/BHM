cmake_minimum_required (VERSION 2.8)

set(CMAKE_MODULE_PATH ..)

find_package(GSL REQUIRED)
message(STATUS "GSL include dir:" ${GSL_INCLUDE_DIR})
message(STATUS "GSL libraries:" ${GSL_LIBRARIES})
message(STATUS "GSL link directories:" ${GSL_LINK_DIRECTORIES})

list(APPEND CMAKE_INSTALL_RPATH ${GSL_LINK_DIRECTORIES})
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

include_directories(
    ${GSL_INCLUDE_DIR}
    . ..
    )
    
add_library (samplinglibrary basic.cpp spline.cpp slot.cpp histogram.cpp)
add_executable(main main.cpp print_spline_grid.cpp)

set_property(TARGET main PROPERTY OUTPUT_NAME "bhm")

add_executable(generator generator.cpp)

target_link_libraries(main samplinglibrary ${GSL_LIBRARIES} libiniparser)
target_link_libraries(generator samplinglibrary ${GSL_LIBRARIES} libiniparser)


install(TARGETS main generator RUNTIME DESTINATION bin)
