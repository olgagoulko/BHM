cmake_minimum_required(VERSION 2.8.12)
project(iniparser C)

add_library(libiniparser src/dictionary.c src/iniparser.c)
target_include_directories(libiniparser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
