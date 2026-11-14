cmake_minimum_required(VERSION 3.16)
project(small_gicp_install_consumer LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(small_gicp CONFIG REQUIRED)
find_package(robot_core CONFIG REQUIRED)
include(${robot_core_DIR}/robot_coreTests.cmake)
enable_testing()
robot_core_add_test(installed_native_helper ../helper_native_test.cpp)
target_link_libraries(installed_native_helper small_gicp::small_gicp)
set_tests_properties(installed_native_helper PROPERTIES ENVIRONMENT
  "OMP_NUM_THREADS=2;OMP_THREAD_LIMIT=2;LD_LIBRARY_PATH=$<TARGET_FILE_DIR:small_gicp::small_gicp>:$ENV{LD_LIBRARY_PATH}")

robot_core_add_test(installed_registration_preparation
  ../registration_preparation_native_test.cpp)
target_link_libraries(installed_registration_preparation small_gicp::small_gicp)
set_tests_properties(installed_registration_preparation PROPERTIES ENVIRONMENT
  "OMP_NUM_THREADS=2;OMP_THREAD_LIMIT=2;LD_LIBRARY_PATH=$<TARGET_FILE_DIR:small_gicp::small_gicp>:$ENV{LD_LIBRARY_PATH}")

robot_core_add_test(installed_loop_registration
  ../loop_registration_native_test.cpp)
target_link_libraries(installed_loop_registration small_gicp::small_gicp)
set_tests_properties(installed_loop_registration PROPERTIES
  RUN_SERIAL TRUE TIMEOUT 90
  ENVIRONMENT "OMP_NUM_THREADS=1;OMP_THREAD_LIMIT=1;OPENBLAS_NUM_THREADS=1;LD_LIBRARY_PATH=$<TARGET_FILE_DIR:small_gicp::small_gicp>:$ENV{LD_LIBRARY_PATH}")

if(UNIX)
  robot_core_add_test(installed_self_fit_tool ../self_fit_tool_native_test.cpp)
  target_link_libraries(installed_self_fit_tool robot_core::alg)
  target_compile_definitions(installed_self_fit_tool PRIVATE
    SELF_FIT_TOOL_PATH="$<TARGET_FILE:small_gicp::sc_self_fit>")
  set_tests_properties(installed_self_fit_tool PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=2;OMP_THREAD_LIMIT=2")
endif()
