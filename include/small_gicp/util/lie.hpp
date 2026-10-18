// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math/se3_exp.h>

namespace small_gicp {

/// @brief Create skew symmetric matrix
/// @param x  Rotation vector
/// @return   Skew symmetric matrix
inline Eigen::Matrix3d skew(const Eigen::Vector3d& x) {
  Eigen::Matrix3d skew = Eigen::Matrix3d::Zero();
  skew(0, 1) = -x[2];
  skew(0, 2) = x[1];
  skew(1, 0) = x[2];
  skew(1, 2) = -x[0];
  skew(2, 0) = -x[1];
  skew(2, 1) = x[0];

  return skew;
}

/*
 * SO3 expmap code taken from Sophus
 * https://github.com/strasdat/Sophus/blob/593db47500ea1a2de5f0e6579c86147991509c59/sophus/so3.hpp#L585
 *
 * Copyright 2011-2017 Hauke Strasdat
 *           2012-2017 Steven Lovegrove
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights  to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

/// @brief SO3 expmap.
/// @param omega  [rx, ry, rz]
/// @return       Quaternion
inline Eigen::Quaterniond so3_exp(const Eigen::Vector3d& omega) {
  const auto q = earth::mantle::math::so3Exp({{omega.x(), omega.y(), omega.z()}});
  return Eigen::Quaterniond(q[3], q[0], q[1], q[2]);
}

// Rotation-first
/// @brief SE3 expmap (Rotation-first).
/// @param a  Twist vector [rx, ry, rz, tx, ty, tz]
/// @return   SE3 matrix
inline Eigen::Isometry3d se3_exp(const Eigen::Matrix<double, 6, 1>& a) {
  const auto pose = earth::mantle::math::se3Exp({{a[0], a[1], a[2], a[3], a[4], a[5]}});
  Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
  for (int row = 0; row < 3; ++row) {
    result.translation()[row] = pose.translation[static_cast<std::size_t>(row)];
    for (int col = 0; col < 3; ++col)
      result.linear()(row, col) = pose.rotation[static_cast<std::size_t>(row * 3 + col)];
  }
  return result;
}

}  // namespace small_gicp
