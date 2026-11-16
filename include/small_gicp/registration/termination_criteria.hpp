#pragma once

#include <small_gicp/registration/linear_system.hpp>
#include <cmath>

namespace small_gicp {

/// @brief Registration termination criteria
struct TerminationCriteria {
  /// @brief Constructor
  TerminationCriteria() : translation_eps(1e-3), rotation_eps(0.1 * M_PI / 180.0) {}

  /// @brief Check the convergence
  /// @param delta  Transformation update vector
  /// @return       True if converged
  bool converged(const factor_detail::Gradient& delta) const { return detail::rotation_step_norm(delta) <= rotation_eps && detail::translation_step_norm(delta) <= translation_eps; }

  double translation_eps;  ///< Translation tolerance [m]
  double rotation_eps;     ///< Rotation tolerance [rad]
};

}  // namespace small_gicp
