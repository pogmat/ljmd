cmake_minimum_required(VERSION 3.10)

if(ENABLE_MPI)


find_package(MPI REQUIRED) 




#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#target library name and sources
add_library(ljmd_mpi_io_lib mpi_io.c)


target_link_libraries(ljmd_mpi_io_lib PRIVATE 
					${MPI_C_LIBRARIES}
					${MPI_C_LINK_FLAGS} 
					ljmd_io_lib
					)
target_compile_definitions(ljmd_mpi_io_lib PRIVATE MPI_ENABLED=1 )

endif()