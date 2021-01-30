cmake_minimum_required(VERSION 3.10)

#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)

#target library name and sources
add_library(ljmd_common_lib common.c)
