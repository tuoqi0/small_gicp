cmake_minimum_required(VERSION 3.16)
project(small_gicp_install_consumer LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(SMALL_GICP_CONSUMER_COMPONENTS "" CACHE STRING "Optional installed components to verify")

# Exercise a native import followed by an explicit optional-component request.
find_package(small_gicp CONFIG REQUIRED)
if(SMALL_GICP_CONSUMER_COMPONENTS)
  find_package(small_gicp CONFIG REQUIRED COMPONENTS ${SMALL_GICP_CONSUMER_COMPONENTS})
endif()
find_package(robot_core CONFIG REQUIRED)
include(${robot_core_DIR}/robot_coreTests.cmake)
enable_testing()

# Reuse the actual library contracts; do not give consumers source include paths.
robot_core_add_test(installed_native_helper ../../src/test/helper_native_test.cpp)
target_link_libraries(installed_native_helper small_gicp::small_gicp)
if("eigen" IN_LIST SMALL_GICP_CONSUMER_COMPONENTS)
  robot_core_add_test(installed_eigen_features ../../src/test/local_features_parity_test.cpp)
  target_link_libraries(installed_eigen_features small_gicp::eigen)
endif()
if("pcl" IN_LIST SMALL_GICP_CONSUMER_COMPONENTS)
  robot_core_add_test(installed_pcl_points ../../src/test/points_test.cpp)
  target_link_libraries(installed_pcl_points small_gicp::pcl)
endif()
if("benchmark" IN_LIST SMALL_GICP_CONSUMER_COMPONENTS)
  add_executable(installed_benchmark_headers benchmark_headers.cpp)
  target_link_libraries(installed_benchmark_headers small_gicp::benchmark)
  add_test(NAME installed_benchmark_headers COMMAND installed_benchmark_headers)
endif()

# A sourced ROS overlay may contain another library with the same SONAME.
get_property(consumer_tests DIRECTORY PROPERTY TESTS)
foreach(consumer_test IN LISTS consumer_tests)
  set_tests_properties(${consumer_test} PROPERTIES ENVIRONMENT
    "OMP_NUM_THREADS=2;LD_LIBRARY_PATH=$<TARGET_FILE_DIR:small_gicp::small_gicp>:$ENV{LD_LIBRARY_PATH}")
endforeach()
