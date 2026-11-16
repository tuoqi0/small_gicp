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
option(BUILD_WITH_MARCH_NATIVE "Build with -march=native" OFF)
include(GNUInstallDirs)

find_package(robot_core CONFIG REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
if(BUILD_WITH_MARCH_NATIVE)
  add_compile_options(-march=native)
endif()
if(MSVC)
  add_compile_definitions(_USE_MATH_DEFINES)
  add_compile_options(/bigobj /openmp:llvm)
endif()

add_library(small_gicp
  src/registration_helper.cpp
  src/registration_preparation.cpp
  src/loop_registration.cpp)
target_include_directories(small_gicp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(small_gicp PUBLIC robot_core::alg OpenMP::OpenMP_CXX)


install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS small_gicp EXPORT small_gicp-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
set(CMAKE_CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/small_gicp"
  CACHE PATH "Install directory for small_gicp CMake package")
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
