# Optional source-tree programs consume PCL through a target, not global paths.
find_package(PCL REQUIRED)
add_library(small_gicp_pcl_dependencies INTERFACE)
target_include_directories(small_gicp_pcl_dependencies SYSTEM INTERFACE ${PCL_INCLUDE_DIRS})
target_compile_options(small_gicp_pcl_dependencies INTERFACE ${PCL_DEFINITIONS})
target_link_libraries(small_gicp_pcl_dependencies INTERFACE ${PCL_LIBRARIES})
