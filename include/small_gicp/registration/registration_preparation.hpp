#pragma once

#include <array>
#include <vector>
#include <small_gicp/points/point_cloud.hpp>

namespace small_gicp {

// XYZ and KNN covariance retain float precision in native double buffers.
// Geometric normals use a physical radius and face the coordinate origin.
// Sampling uses the OMP sampler's default workers; threads controls covariance
// and radius-normal estimation. Empty input returns a nonnull empty cloud.
PointCloud::Ptr prepare_registration_cloud(
  const std::vector<std::array<float, 3>>& coordinates, double leaf,
  int neighbors, int threads, double normal_radius);

}  // namespace small_gicp
