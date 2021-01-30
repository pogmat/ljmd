cmake_minimum_required(VERSION 3.10)

include_directories(${ljmd_SOURCE_DIR}/include)

set(ljmd_SOURCE_FILES ${ljmd_SOURCE_DIR}/src )

message("-- Compiling unit tests")



set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(GTest)

add_executable(test_common ${ljmd_SOURCE_FILES}/common/test_common.cpp)
#add_executable(test_common common/test_common.cpp)
target_include_directories(test_common PRIVATE ${CMAKE_MODULE_PATH}/src)
target_link_libraries(test_common PRIVATE GTest::GTestMain GTest::GTest ljmd_common_lib)
add_test(NAME common COMMAND test_common)
  
