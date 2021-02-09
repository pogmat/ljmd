cmake_minimum_required(VERSION 3.10)

#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)

#target library name and sources
add_library(ljmd_phys_lib integration.c force.c ekin.c utils.c cell.c)
target_link_libraries(ljmd_phys_lib PRIVATE m)

#list of compilation options
list(APPEND phys_flags -Wall -g -O3 -ffast-math -fomit-frame-pointer)

#compile options
target_compile_options (ljmd_phys_lib PRIVATE ${phys_flags}  )

get_target_property(PHYS_COMPILE_OPT ljmd_phys_lib COMPILE_OPTIONS)
message("-- ljmd_phys_lib compiler flags are: ${PHYS_COMPILE_OPT}")
