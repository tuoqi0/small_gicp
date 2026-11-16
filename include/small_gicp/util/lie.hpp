#pragma once

#include <math/fixed_matrix_ops.h>
#include <math/se3_exp.h>

namespace small_gicp {

// Homogeneous source-to-target transform. Storage and products use Mantle.
using Transform = earth::mantle::collection::Matrix<double, 4, 4>;

/// @brief SE3 exponential, rotation-first increment [rx, ry, rz, tx, ty, tz].
inline Transform se3_exp(const earth::mantle::collection::Matrix<double, 6, 1>& delta) {
  const auto pose = earth::mantle::math::se3Exp({{delta[0][0], delta[1][0], delta[2][0],
      delta[3][0], delta[4][0], delta[5][0]}});
  auto transform = earth::mantle::math::identityMatrix<double, 4>();
  for (size_t row = 0; row < 3; ++row) {
    for (size_t col = 0; col < 3; ++col) transform[row][col] = pose.rotation[row * 3 + col];
    transform[row][3] = pose.translation[row];
  }
  return transform;
}

inline earth::mantle::collection::Matrix<double, 3, 3> rotation_matrix(const Transform& transform) {
  earth::mantle::collection::Matrix<double, 3, 3> rotation{};
  for (size_t row = 0; row < 3; ++row)
    for (size_t col = 0; col < 3; ++col) rotation[row][col] = transform[row][col];
  return rotation;
}

}  // namespace small_gicp
