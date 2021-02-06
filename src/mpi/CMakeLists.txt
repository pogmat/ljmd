cmake_minimum_required(VERSION 3.10)

find_package(MPI REQUIRED) 

message("fount MPI C compiler : ${MPI_C_COMPILER}")
message("fount MPI C libs : ${MPI_C_LIBRARIES}")



#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#target library name and sources
add_library(ljmd_mpi_lib mpi_utils.c )


target_link_libraries(ljmd_mpi_lib PRIVATE 
					${MPI_C_LIBRARIES}
					${MPI_C_LINK_FLAGS} 
					)

