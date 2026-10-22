// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <math/fixed_matrix_ops.h>
#include <math/se3_exp.h>

namespace small_gicp {

// Homogeneous source-to-target transform. Storage and products use Mantle.
using Transform = earth::mantle::collection::Matrix<double, 4, 4>;

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
