cmake_minimum_required(VERSION 3.10)

find_package(MPI REQUIRED) 


#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#target library name and sources
add_library(ljmd_mpi_comm_lib mpi_comm.c )


target_link_libraries(ljmd_mpi_comm_lib PRIVATE 
					${MPI_C_LIBRARIES}
					${MPI_C_LINK_FLAGS} 
					)

