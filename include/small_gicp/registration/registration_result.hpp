#pragma once

#include <small_gicp/factors/linearization.hpp>
#include <small_gicp/util/lie.hpp>

namespace small_gicp {

/// @brief Registration result
struct RegistrationResult {
  RegistrationResult(const Transform& T = earth::mantle::math::identityMatrix<double, 4>())
  : T_target_source(T),
    converged(false),
    iterations(0),
    num_inliers(0),
    H{},
    b{},
    error(0.0) {}

  Transform T_target_source;  ///<  Estimated transformation

  bool converged;      ///< If the optimization converged
  size_t iterations;   ///< Number of optimization iterations
  size_t num_inliers;  ///< Number of inliear points

  factor_detail::Hessian H;  ///< Final information matrix
  factor_detail::Gradient b;  ///< Final information vector
  double error;                   ///< Final error
};

}  // namespace small_gicp
