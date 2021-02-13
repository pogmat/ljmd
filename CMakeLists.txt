cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(ljmd VERSION 1.1 LANGUAGES C)

#specify target exe name and source
add_executable(serial  src/ljmd.c)
set_target_properties(serial PROPERTIES OUTPUT_NAME "ljmd-serial.x")

#subdirectories to run their own cmake files 
list(APPEND sub_dirs 
					${ljmd_SOURCE_DIR}/src/common  
					${ljmd_SOURCE_DIR}/src/io  
					${ljmd_SOURCE_DIR}/src/physics 
					)

if(ENABLE_TESTING)
					enable_testing()
					enable_language(CXX)
					list(APPEND sub_dirs ${ljmd_SOURCE_DIR}/tests )
endif()

foreach(subdir ${sub_dirs})
  					add_subdirectory(${subdir})
endforeach()


#include dir for the header files
include_directories(${ljmd_SOURCE_DIR}/include)



#link libraries and the lib source files
target_link_libraries(serial PRIVATE 
					ljmd_phys_lib
					ljmd_io_lib
					ljmd_common_lib
					)





#compilation options list
list(APPEND main_flags -Wall -g -O3 -ffast-math -fomit-frame-pointer)

#compile options
target_compile_options (serial PRIVATE ${main_flags}  )
target_compile_definitions(serial PRIVATE LJMD_VERSION=1.0)

get_target_property(LJMD_COMPILE_OPT serial COMPILE_OPTIONS)
message("-- ljmd main compiler flags are: ${LJMD_COMPILE_OPT}")



##########################
# _OMP_NAIVE
##########################

#specify target exe name and source
add_executable(omp_naive src/ljmd.c)
set_target_properties(omp_naive PROPERTIES OUTPUT_NAME "ljmd_omp_naive.x")

#link libraries and the lib source files
find_package(OpenMP REQUIRED)
target_link_libraries(omp_naive PRIVATE ljmd_phys_lib_omp_naive	ljmd_io_lib ljmd_common_lib OpenMP::OpenMP_C)

#compile options
target_compile_options (omp_naive PRIVATE ${main_flags})
target_compile_definitions(omp_naive PRIVATE LJMD_VERSION=1.0 _OMP_NAIVE)