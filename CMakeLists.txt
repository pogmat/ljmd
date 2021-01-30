cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(ljmd VERSION 1.1 LANGUAGES C)

#specify target exe name and source
add_executable(serial  src/ljmd.c)
set_target_properties(serial PROPERTIES OUTPUT_NAME "ljmd-serial.x")

#subdirectories to run their own cmake files 
add_subdirectory(${ljmd_SOURCE_DIR}/src/common)
add_subdirectory(${ljmd_SOURCE_DIR}/src/io)
add_subdirectory(${ljmd_SOURCE_DIR}/src/physics)
add_subdirectory(${ljmd_SOURCE_DIR}/tests)


#include dir for the header files
include_directories(${ljmd_SOURCE_DIR}/include)



#link libraries and the lib source files
target_link_libraries(serial PRIVATE 
					ljmd_phys_lib
					ljmd_io_lib
					ljmd_common_lib
					)

link_directories(
				${ljmd_SOURCE_DIR}/src/common
				${ljmd_SOURCE_DIR}/src/io
				${ljmd_SOURCE_DIR}/src/physics
				)



#compilation options list
list(APPEND main_flags -Wall -g -O3 -ffast-math -fomit-frame-pointer)

#compile options
target_compile_options (serial PRIVATE ${main_flags}  )
target_compile_definitions(serial PRIVATE LJMD_VERSION=1.0)

get_target_property(LJMD_COMPILE_OPT serial COMPILE_OPTIONS)
message("-- ljmd main compiler flags are: ${LJMD_COMPILE_OPT}")





