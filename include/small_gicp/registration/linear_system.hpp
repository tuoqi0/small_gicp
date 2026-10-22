// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <small_gicp/factors/linearization.hpp>

namespace small_gicp { namespace detail {

// Eigen's original LDLT reads the lower triangle. Accumulation can leave the
// two triangles slightly different; preserve that convention and -gradient.
inline bool solve_damped_system(const factor_detail::Hessian& hessian,
                                const factor_detail::Gradient& gradient,
                                double damping, factor_detail::Gradient& delta) {
  std::array<std::array<double, 6>, 6> matrix{};
  std::array<double, 6> rhs{}, solution{};
  for (std::size_t row = 0; row < 6; ++row) {
    rhs[row] = -gradient[row][0];
    for (std::size_t col = 0; col < 6; ++col)
      matrix[row][col] = hessian[std::max(row, col)][std::min(row, col)];
    matrix[row][row] += damping;
  }
  if (!earth::mantle::math::solvePositiveDefinite(matrix, rhs, solution)) return false;
  delta = factor_detail::Gradient::from(solution.data(), solution.size());
  return true;
}

// Twist order is rotation first, followed by translation.
inline double rotation_step_norm(const factor_detail::Gradient& delta) {
  return std::sqrt(delta[0][0] * delta[0][0] + delta[1][0] * delta[1][0] + delta[2][0] * delta[2][0]);
}
inline double translation_step_norm(const factor_detail::Gradient& delta) {
  return std::sqrt(delta[3][0] * delta[3][0] + delta[4][0] * delta[4][0] + delta[5][0] * delta[5][0]);
}

}}  // namespace small_gicp::detail
