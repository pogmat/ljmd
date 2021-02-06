cmake_minimum_required(VERSION 3.10)

include_directories(${ljmd_SOURCE_DIR}/include)

set(ljmd_SOURCE_FILES ${ljmd_SOURCE_DIR}/src )

message("-- Compiling unit tests")



set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(GTest)

#common.c tests

add_executable(test_common ${ljmd_SOURCE_FILES}/common/test_common.cpp)
set_target_properties(test_common PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
target_include_directories(test_common PRIVATE ${CMAKE_MODULE_PATH}/src)
target_link_libraries(test_common PRIVATE GTest::GTestMain GTest::GTest ljmd_common_lib)
add_test(NAME common COMMAND test_common )
  
# io tests
add_executable(test_input ${ljmd_SOURCE_FILES}/io/test_input.cpp)
set_target_properties(test_input PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
target_include_directories(test_input PRIVATE ${CMAKE_MODULE_PATH}/src ${ljmd_SOURCE_DIR}/src/io ${ljmd_SOURCE_DIR}/src/common)
target_link_libraries(test_input PRIVATE GTest::GTestMain GTest::GTest ljmd_io_lib)
add_test(NAME input COMMAND test_input)


# utils tests

add_executable(test_utils ${ljmd_SOURCE_FILES}/physics/test_utils.cpp)
set_target_properties(test_utils PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
target_include_directories(test_utils PRIVATE ${CMAKE_MODULE_PATH}/src ${ljmd_SOURCE_DIR}/src/physics)
target_link_libraries(test_utils PRIVATE GTest::GTestMain GTest::GTest ljmd_phys_lib)
add_test(NAME utils COMMAND test_utils )



# force tests

add_executable(test_force ${ljmd_SOURCE_FILES}/physics/test_force.cpp)
set_target_properties(test_force PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
target_include_directories(test_force PRIVATE ${CMAKE_MODULE_PATH}/src ${ljmd_SOURCE_DIR}/src/physics)
target_link_libraries(test_force PRIVATE GTest::GTestMain GTest::GTest ljmd_phys_lib ljmd_common_lib)
add_test(NAME force COMMAND test_force )


# ekin and temp tests

add_executable(test_ek_t ${ljmd_SOURCE_FILES}/physics/test_ekin.cpp)
set_target_properties(test_ek_t PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
target_include_directories(test_ek_t PRIVATE ${CMAKE_MODULE_PATH}/src ${ljmd_SOURCE_DIR}/src/physics)
target_link_libraries(test_ek_t PRIVATE GTest::GTestMain GTest::GTest ljmd_phys_lib ljmd_common_lib)
add_test(NAME ekin COMMAND test_ek_t )

# verlet integration tests

add_executable(test_verlet ${ljmd_SOURCE_FILES}/physics/test_integration.cpp)
set_target_properties(test_verlet PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
target_include_directories(test_verlet PRIVATE ${CMAKE_MODULE_PATH}/src ${ljmd_SOURCE_DIR}/src/physics)
target_link_libraries(test_verlet PRIVATE GTest::GTestMain GTest::GTest ljmd_phys_lib ljmd_common_lib)
add_test(NAME verlet COMMAND test_verlet )

# reference test
find_program(BASH_PROGRAM bash)
find_program(AWK_PROGRAM awk)
find_program(CMP_PROGRAM cmp)
if(BASH_PROGRAM AND AWK_PROGRAM AND CMP_PROGRAM)
    add_test(NAME comparison COMMAND ${ljmd_SOURCE_DIR}/tests/compare.sh $<TARGET_FILE:serial> ${CMAKE_MODULE_PATH})
endif()





#mpi tests

if(ENABLE_MPI)
	find_package(MPI REQUIRED) 
	

	message("fount MPI C compiler : ${MPI_C_COMPILER}")
	message("fount MPI CXX compiler : ${MPI_CXX_COMPILER}")
	
	function(create_test_mpi)
		
		add_executable(test_mpi ${ljmd_SOURCE_FILES}/mpi/test_mpi.cpp)
		set_target_properties(test_mpi PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${ljmd_SOURCE_DIR}/build)
		target_include_directories(test_mpi PRIVATE ${MPI_CXX_INCLUDE_PATH} )
		target_link_libraries(test_mpi PRIVATE GTest::GTestMain GTest::GTest ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS}  )
		
		set(TEST_MPI_NPROCS 4)
		
		cmake_parse_arguments(TEST_MPI_ "" "NPROCS" "COMMAND" ${ARGN})
		message("${TEST_MPI_NAME} command : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_MPI_NPROCS} ${MPIEXEC_PREFLAGS} ./test_mpi ${TEST_MPI__COMMAND}")
		
		add_test(
			NAME mpi 
			COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_MPI_NPROCS} ${TEST_MPI__COMMAND} ${TEST_MPI_NPROCS}
		)
		
		
	
	endfunction(create_test_mpi)
	
	create_test_mpi(NPROCS 4 COMMAND $<TARGET_FILE:test_mpi> )
	
endif()

