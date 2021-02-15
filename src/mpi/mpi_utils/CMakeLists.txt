cmake_minimum_required(VERSION 3.10)

if(ENABLE_MPI)

#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#target library name and sources
add_library(ljmd_mpi_util_lib mpi_utils.c )


target_link_libraries(ljmd_mpi_util_lib PRIVATE m )
target_compile_definitions(ljmd_mpi_util_lib PRIVATE MPI_ENABLED=1 )

endif()