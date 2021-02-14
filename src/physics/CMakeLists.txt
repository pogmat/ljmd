cmake_minimum_required(VERSION 3.10)

#header files directory
include_directories(${ljmd_SOURCE_DIR}/include)

#target library name and sources
add_library(ljmd_phys_lib integration.c force.c ekin.c utils.c)
target_link_libraries(ljmd_phys_lib PRIVATE m)

#list of compilation options
list(APPEND phys_flags -Wall -g -O3 -ffast-math -fomit-frame-pointer)

#compile options
target_compile_options (ljmd_phys_lib PRIVATE ${phys_flags}  )

get_target_property(PHYS_COMPILE_OPT ljmd_phys_lib COMPILE_OPTIONS)
message("-- ljmd_phys_lib compiler flags are: ${PHYS_COMPILE_OPT}")

##########################
# _OMP_NAIVE
##########################

# target library name and sources (except force.c)
add_library(ljmd_phys_lib_omp_naive integration.c ekin.c utils.c)
find_package(OpenMP REQUIRED)
target_link_libraries(ljmd_phys_lib_omp_naive PRIVATE m OpenMP::OpenMP_C)

# compile options (except force.c)
target_compile_options(ljmd_phys_lib_omp_naive PRIVATE ${phys_flags})
target_compile_definitions(ljmd_phys_lib_omp_naive PRIVATE _OMP_NAIVE)

# target library name and sources (for force.c)
add_library(ljmd_force_lib_omp_naive force.c)
target_link_libraries(ljmd_force_lib_omp_naive PRIVATE m OpenMP::OpenMP_C)

# compile options (for force.c)
target_compile_options(ljmd_force_lib_omp_naive PRIVATE ${phys_flags})
target_compile_definitions(ljmd_force_lib_omp_naive PRIVATE _OMP_NAIVE)

###########################
# _OMP_3RD_LAW
###########################

# target library name and sources (for force.c)
add_library(ljmd_force_lib_omp_3rd_law force.c)
target_link_libraries(ljmd_force_lib_omp_3rd_law PRIVATE m OpenMP::OpenMP_C)

# compile options (for force.c)
target_compile_options(ljmd_force_lib_omp_3rd_law PRIVATE ${phys_flags})
target_compile_definitions(ljmd_force_lib_omp_3rd_law PRIVATE _OMP_NAIVE _OMP_3RD_LAW)
