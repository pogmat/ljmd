cmake_minimum_required(VERSION 3.10)

#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)

if(ENABLE_MPI)
	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif()

#target library name and sources
add_library(ljmd_io_lib input.c output.c)

if(ENABLE_MPI)
	target_link_libraries(ljmd_io_lib PRIVATE 
					${MPI_C_LIBRARIES}
					${MPI_C_LINK_FLAGS} 
					ljmd_mpi_util_lib
					ljmd_mpi_comm_lib
					)
	target_compile_definitions(ljmd_io_lib PRIVATE MPI_ENABLED=1 )
endif()

