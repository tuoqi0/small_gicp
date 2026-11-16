#pragma once

#include <array>
#include <cstddef>
#include <collection/fixed_matrix.h>
#include <math/fixed_matrix_ops.h>
#include <math/positive_definite_solve.h>

namespace small_gicp { namespace factor_detail {

using Point = earth::mantle::collection::Matrix<double, 4, 1>;
using Rotation = earth::mantle::collection::Matrix<double, 3, 3>;
using Covariance = earth::mantle::collection::Matrix<double, 4, 4>;
using Jacobian = earth::mantle::collection::Matrix<double, 4, 6>;
using Hessian = earth::mantle::collection::Matrix<double, 6, 6>;
using Gradient = earth::mantle::collection::Matrix<double, 6, 1>;

// Rotation-first perturbation. The homogeneous coordinate has no derivative.
inline Jacobian pointJacobian(const Point& source, const Rotation& rotation) {
  using earth::mantle::math::operator*;
  const auto xyz = earth::mantle::collection::Matrix<double, 3, 1>::from(
      {source[0][0], source[1][0], source[2][0]});
  const auto orientation = rotation * earth::mantle::math::hat(xyz);
  Jacobian result{};
  for (std::size_t row = 0; row < 3; ++row) {
    for (std::size_t col = 0; col < 3; ++col) {
      result[row][col] = orientation[row][col];
      result[row][col + 3] = -rotation[row][col];
    }
  }
  return result;
}

inline double squaredError(const Point& residual) {
  using earth::mantle::math::operator*;
  return (0.5 * earth::mantle::math::transpose(residual) * residual)[0][0];
}

inline double squaredError(const Point& residual, const Covariance& precision) {
  using earth::mantle::math::operator*;
  return (0.5 * earth::mantle::math::transpose(residual) * precision * residual)[0][0];
}

inline void linearizeResidual(const Jacobian& jacobian, const Point& residual,
    Hessian& hessian, Gradient& gradient, double& error) {
  using earth::mantle::math::operator*;
  const auto transposed = earth::mantle::math::transpose(jacobian);
  hessian = transposed * jacobian;
  gradient = transposed * residual;
  error = squaredError(residual);
}

inline void linearizeResidual(const Jacobian& jacobian, const Point& residual,
    const Covariance& precision, Hessian& hessian, Gradient& gradient, double& error) {
  using earth::mantle::math::operator*;
  const auto weighted_transpose = earth::mantle::math::transpose(jacobian) * precision;
  hessian = weighted_transpose * jacobian;
  gradient = weighted_transpose * residual;
  error = squaredError(residual, precision);
}

// Keep the existing SPD requirement. A failed solve must not publish a partial weight.
inline bool gicpPrecision(const Covariance& target, const Covariance& source,
    const Covariance& transform, Covariance& output) {
  using earth::mantle::math::operator*;
  using earth::mantle::math::operator+;
  const auto combined = target + transform * source * earth::mantle::math::transpose(transform);
  std::array<std::array<double, 3>, 3> covariance{}, identity{}, precision{};
  for (std::size_t row = 0; row < 3; ++row) {
    identity[row][row] = 1;
    for (std::size_t col = 0; col < 3; ++col) covariance[row][col] = combined[row][col];
  }
  if (!earth::mantle::math::solvePositiveDefinite(covariance, identity, precision)) return false;
  Covariance weight{};
  for (std::size_t row = 0; row < 3; ++row)
    for (std::size_t col = 0; col < 3; ++col) weight[row][col] = precision[row][col];
  output = weight;
  return true;
}

}}  // namespace small_gicp::factor_detail
