// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <small_gicp/factors/linearization.hpp>
#include <Eigen/Geometry>

namespace small_gicp {

/// @brief Registration result
struct RegistrationResult {
  RegistrationResult(const Eigen::Isometry3d& T = Eigen::Isometry3d::Identity())
  : T_target_source(T),
    converged(false),
    iterations(0),
    num_inliers(0),
    H{},
    b{},
    error(0.0) {}

  Eigen::Isometry3d T_target_source;  ///<  Estimated transformation

  bool converged;      ///< If the optimization converged
  size_t iterations;   ///< Number of optimization iterations
  size_t num_inliers;  ///< Number of inliear points

  factor_detail::Hessian H;  ///< Final information matrix
  factor_detail::Gradient b;  ///< Final information vector
  double error;                   ///< Final error
};

}  // namespace small_gicp
