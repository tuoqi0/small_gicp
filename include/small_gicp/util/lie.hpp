// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math/se3_exp.h>
#include <collection/fixed_matrix.h>

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
  const auto q = earth::mantle::math::so3Exp({{omega[0], omega[1], omega[2]}});
  // Mantle stores x/y/z/w; Eigen's constructor takes w/x/y/z.
  return Eigen::Quaterniond(q[3], q[0], q[1], q[2]);
}

// Rotation-first
/// @brief SE3 expmap (Rotation-first).
/// @param a  Twist vector [rx, ry, rz, tx, ty, tz]
/// @return   SE3 matrix
inline Eigen::Isometry3d se3_exp(const Eigen::Matrix<double, 6, 1>& a) {
  const auto native = earth::mantle::math::se3Exp({{a[0], a[1], a[2], a[3], a[4], a[5]}});
  Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
  result.linear() = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(
      native.rotation.data());
  result.translation() = Eigen::Map<const Eigen::Vector3d>(native.translation.data());
  return result;
}

// Temporary Eigen pose boundary; the optimizer increment remains native.
inline Eigen::Isometry3d se3_exp(const earth::mantle::collection::Matrix<double, 6, 1>& a) {
  const auto native = earth::mantle::math::se3Exp({{a[0][0], a[1][0], a[2][0], a[3][0], a[4][0], a[5][0]}});
  Eigen::Isometry3d result = Eigen::Isometry3d::Identity();
  result.linear() = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(native.rotation.data());
  result.translation() = Eigen::Map<const Eigen::Vector3d>(native.translation.data());
  return result;
}

}  // namespace small_gicp
