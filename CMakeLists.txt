cmake_minimum_required(VERSION 3.16)
project(small_gicp VERSION 1.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug MinSizeRel RelWithDebInfo)
endif()

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(CMAKE_POSITION_INDEPENDENT_CODE "Generate position-independent code" ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
option(BUILD_HELPER "Build the native registration helper library" ON)
option(BUILD_WITH_TBB "Enable native TBB headers and their tests" OFF)
option(BUILD_WITH_MARCH_NATIVE "Build with -march=native" OFF)
option(ENABLE_COVERAGE "Enable native test coverage" OFF)
set(BUILD_WITH_OPENMP auto CACHE STRING "Build with OpenMP")
include(CTest)

find_package(robot_core CONFIG REQUIRED)
if(BUILD_WITH_OPENMP STREQUAL "auto")
  find_package(OpenMP)
  set(BUILD_WITH_OPENMP ${OpenMP_FOUND})
elseif(BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
endif()
if(BUILD_WITH_TBB)
  find_package(TBB REQUIRED)
  add_compile_definitions(BUILD_WITH_TBB)
endif()
if(BUILD_WITH_MARCH_NATIVE)
  add_compile_options(-march=native)
endif()
if(MSVC)
  add_compile_definitions(_USE_MATH_DEFINES)
  add_compile_options(/bigobj)
  if(BUILD_WITH_OPENMP)
    add_compile_options(/openmp:llvm)
  endif()
endif()

if(ENABLE_COVERAGE)
  set(CMAKE_CXX_FLAGS "-O0 -coverage")
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml REQUIRED)
  add_custom_target(coverage
    COMMAND ${LCOV} --directory . --capture --output-file coverage.info --ignore-errors mismatch
    COMMAND ${LCOV} --remove coverage.info -o coverage.info '/usr/*' --ignore-errors mismatch
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if(BUILD_TESTING)
  set(BUILD_HELPER ON CACHE BOOL "Native contracts require the helper library" FORCE)
endif()
if(BUILD_HELPER)
  add_library(small_gicp
    src/small_gicp/registration/registration.cpp
    src/small_gicp/registration/registration_helper.cpp)
  target_include_directories(small_gicp PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
  target_link_libraries(small_gicp PUBLIC
    robot_core::alg $<TARGET_NAME_IF_EXISTS:OpenMP::OpenMP_CXX>)
endif()

if(BUILD_TESTING)
  find_package(robot_core CONFIG REQUIRED)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  include(${robot_core_DIR}/robot_coreTests.cmake)
  set(NATIVE_TEST_SOURCES
    src/test/helper_native_test.cpp
    src/test/linear_system_native_test.cpp
    src/test/linearization_native_test.cpp
    src/test/projection_native_test.cpp
    src/test/registration_native_test.cpp
    src/test/point_storage_native_test.cpp
    src/test/flat_container_native_test.cpp
    src/test/ros_storage_native_test.cpp
    src/test/optimizer_system_test.cpp
    src/test/gaussian_voxel_native_test.cpp
    src/test/sort_omp_test.cpp
    src/test/vector_test.cpp)
  if(BUILD_WITH_TBB)
    list(APPEND NATIVE_TEST_SOURCES src/test/sort_tbb_test.cpp)
  endif()
  foreach(TEST_SOURCE IN LISTS NATIVE_TEST_SOURCES)
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    robot_core_add_test(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} small_gicp OpenMP::OpenMP_CXX)
    if(TEST_NAME STREQUAL "sort_tbb_test")
      target_link_libraries(${TEST_NAME} TBB::tbb)
    endif()
    set_tests_properties(${TEST_NAME} PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=2;OMP_THREAD_LIMIT=2;LD_LIBRARY_PATH=$<TARGET_FILE_DIR:small_gicp>:$ENV{LD_LIBRARY_PATH}")
  endforeach()
endif()

include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
if(BUILD_HELPER)
  install(TARGETS small_gicp EXPORT small_gicp-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  set(CMAKE_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/small_gicp"
    CACHE PATH "Install directory for CMake config files")
  include(CMakePackageConfigHelpers)
  install(EXPORT small_gicp-targets FILE small_gicp-targets.cmake
    NAMESPACE small_gicp:: DESTINATION ${CMAKE_CONFIG_INSTALL_DIR})
  configure_package_config_file(cmake/small_gicp-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/small_gicp-config.cmake
    INSTALL_DESTINATION ${CMAKE_CONFIG_INSTALL_DIR})
  write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/small_gicp-config-version.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/small_gicp-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/small_gicp-config-version.cmake
    DESTINATION ${CMAKE_CONFIG_INSTALL_DIR})
endif()
