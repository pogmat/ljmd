cmake_minimum_required(VERSION 3.10)

include_directories(${ljmd_SOURCE_DIR}/include)

if(ENABLE_TESTING)
  enable_testing()
  enable_language(CXX)

  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
  include(GTest)
  
  add_executable(test_common common/test_common.cpp)
  target_include_directories(test_common PRIVATE ${CMAKE_MODULE_PATH}/src)
  target_link_libraries(test_common PRIVATE GTest::GTestMain GTest::GTest ljmd_common_lib)
  add_test(NAME common COMMAND test_common)
  
endif()
