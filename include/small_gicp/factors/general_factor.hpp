// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <small_gicp/util/lie.hpp>
#include <array>
#include <cmath>
#include <small_gicp/factors/linearization.hpp>

namespace small_gicp {

/// @brief Null factor that gives no constraints.
struct NullFactor {
  NullFactor() = default;

  /// @brief Update linearized system consisting of linearized per-point factors.
  /// @param target       Target point cloud
  /// @param source       Source point cloud
  /// @param target_tree  Nearest neighbor search for the target point cloud
  /// @param T            Linearization point
  /// @param H            [in/out] Linearized information matrix.
  /// @param b            [in/out] Linearized information vector.
  /// @param e            [in/out] Error at the linearization point.
  template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree>
  void update_linearized_system(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const Transform& T,
    factor_detail::Hessian* H,
    factor_detail::Gradient* b,
    double* e) const {}

  /// @brief Update error consisting of per-point factors.
  /// @param target   Target point cloud
  /// @param source   Source point cloud
  /// @param T        Evaluation point
  /// @param e        [in/out] Error at the evaluation point.
  template <typename TargetPointCloud, typename SourcePointCloud>
  void update_error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T, double* e) const {}
};

/// @brief Factor to restrict the degrees of freedom of optimization (e.g., fixing roll, pitch rotation).
/// @note  This factor only enables *soft* constraints. The source point cloud can move to the restricted directions slightly).
struct RestrictDoFFactor {
  /// @brief Constructor.
  RestrictDoFFactor() {
    lambda = 1e9;
    mask.fill(1.0);
  }

  /// @brief Set rotation mask. (1.0 = active, 0.0 = inactive)
  template <typename Mask>
  void set_rotation_mask(const Mask& rot_mask) {
    for (size_t axis = 0; axis < 3; ++axis) mask[axis] = rot_mask[axis];
  }

  /// @brief  Set translation mask. (1.0 = active, 0.0 = inactive)
  template <typename Mask>
  void set_translation_mask(const Mask& trans_mask) {
    for (size_t axis = 0; axis < 3; ++axis) mask[axis + 3] = trans_mask[axis];
  }

  /// @brief Update linearized system consisting of linearized per-point factors.
  template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree>
  void update_linearized_system(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const Transform& T,
    factor_detail::Hessian* H,
    factor_detail::Gradient* b,
    double* e) const {
    for (size_t axis = 0; axis < mask.size(); ++axis)
      (*H)[axis][axis] += lambda * std::abs(mask[axis] - 1.0);
  }

  /// @brief Update error consisting of per-point factors.
  template <typename TargetPointCloud, typename SourcePointCloud>
  void update_error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T, double* e) const {}

  double lambda;                    ///< Regularization parameter (Increasing this makes the constraint stronger)
  std::array<double, 6> mask;  ///< Mask for restricting DoF (rx, ry, rz, tx, ty, tz)
};

}  // namespace small_gicp
