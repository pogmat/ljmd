cmake_minimum_required(VERSION 3.10)

if(ENABLE_MPI)

#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)

set(phys_src_path ${ljmd_SOURCE_DIR}/src/physics)

#target library name and sources
add_library(ljmd_mpi_phys_lib 
			mpi_integration.c 
			mpi_force.c 
			mpi_ekin.c 
			${phys_src_path}/utils.c
			)
			
target_link_libraries(ljmd_mpi_phys_lib PRIVATE m)

#list of compilation options
list(APPEND mpi_phys_flags -Wall -g -O3 -ffast-math -fomit-frame-pointer)

#compile options
target_compile_options (ljmd_mpi_phys_lib PRIVATE ${mpi_phys_flags}  )
target_compile_definitions(ljmd_mpi_phys_lib PRIVATE MPI_ENABLED=1 )

get_target_property(PHYS_COMPILE_OPT ljmd_mpi_phys_lib COMPILE_OPTIONS)
message("-- ljmd_mpi_phys_lib compiler flags are: ${PHYS_COMPILE_OPT}")

endif()