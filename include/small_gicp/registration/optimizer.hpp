// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>
#include <algorithm>
#include <small_gicp/registration/linear_system.hpp>
#include <small_gicp/util/lie.hpp>
#include <small_gicp/registration/registration_result.hpp>

namespace small_gicp {

/// @brief GaussNewton optimizer
struct GaussNewtonOptimizer {
  GaussNewtonOptimizer() : verbose(false), max_iterations(20), lambda(1e-6) {}

  template <
    typename TargetPointCloud,
    typename SourcePointCloud,
    typename TargetTree,
    typename CorrespondenceRejector,
    typename TerminationCriteria,
    typename Reduction,
    typename Factor,
    typename GeneralFactor>
  RegistrationResult optimize(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const CorrespondenceRejector& rejector,
    const TerminationCriteria& criteria,
    Reduction& reduction,
    const Transform& init_T,
    std::vector<Factor>& factors,
    GeneralFactor& general_factor) const {
    //
    if (verbose) {
      std::cout << "--- GN optimization ---" << std::endl;
    }

    using earth::mantle::math::operator*;
    RegistrationResult result(init_T);
    for (int i = 0; i < max_iterations && !result.converged; i++) {
      // Linearize
      auto [H, b, e] = reduction.linearize(target, source, target_tree, rejector, result.T_target_source, factors);
      general_factor.update_linearized_system(target, source, target_tree, result.T_target_source, &H, &b, &e);

      // Solve linear system
      factor_detail::Gradient delta{};
      if (!detail::solve_damped_system(H, b, lambda, delta)) {
        result.iterations = i;
        result.H = H;
        result.b = b;
        result.error = e;
        break;
      }

      if (verbose) {
        std::cout << "iter=" << i << " e=" << e << " lambda=" << lambda << " dt=" << detail::translation_step_norm(delta) << " dr=" << detail::rotation_step_norm(delta) << std::endl;
      }

      result.converged = criteria.converged(delta);
      result.T_target_source = result.T_target_source * se3_exp(delta);
      result.iterations = i;
      result.H = H;
      result.b = b;
      result.error = e;
    }

    result.num_inliers = std::count_if(factors.begin(), factors.end(), [](const auto& factor) { return factor.inlier(); });

    return result;
  }

  bool verbose;        ///< If true, print debug messages
  int max_iterations;  ///< Max number of optimization iterations
  double lambda;       ///< Damping factor (Increasing this makes optimization slow but stable)
};

/// @brief LevenbergMarquardt optimizer
struct LevenbergMarquardtOptimizer {
  LevenbergMarquardtOptimizer() : verbose(false), max_iterations(20), max_inner_iterations(10), init_lambda(1e-3), lambda_factor(10.0) {}

  template <
    typename TargetPointCloud,
    typename SourcePointCloud,
    typename TargetTree,
    typename CorrespondenceRejector,
    typename TerminationCriteria,
    typename Reduction,
    typename Factor,
    typename GeneralFactor>
  RegistrationResult optimize(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const CorrespondenceRejector& rejector,
    const TerminationCriteria& criteria,
    Reduction& reduction,
    const Transform& init_T,
    std::vector<Factor>& factors,
    GeneralFactor& general_factor) const {
    //
    if (verbose) {
      std::cout << "--- LM optimization ---" << std::endl;
    }

    double lambda = init_lambda;
    using earth::mantle::math::operator*;
    RegistrationResult result(init_T);
    for (int i = 0; i < max_iterations && !result.converged; i++) {
      // Linearize
      auto [H, b, e] = reduction.linearize(target, source, target_tree, rejector, result.T_target_source, factors);
      general_factor.update_linearized_system(target, source, target_tree, result.T_target_source, &H, &b, &e);

      // Lambda iteration
      bool success = false;
      for (int j = 0; j < max_inner_iterations; j++) {
        // Solve with damping
        factor_detail::Gradient delta{};
        if (!detail::solve_damped_system(H, b, lambda, delta)) {
          lambda *= lambda_factor;
          continue;
        }

        // Validate new solution
        const Transform new_T = result.T_target_source * se3_exp(delta);
        double new_e = reduction.error(target, source, new_T, factors);
        general_factor.update_error(target, source, new_T, &e);

        if (verbose) {
          std::cout << "iter=" << i << " inner=" << j << " e=" << e << " new_e=" << new_e << " lambda=" << lambda << " dt=" << detail::translation_step_norm(delta)
                    << " dr=" << detail::rotation_step_norm(delta) << std::endl;
        }

        if (new_e <= e) {
          // Error decreased, decrease lambda
          result.converged = criteria.converged(delta);
          result.T_target_source = new_T;
          lambda /= lambda_factor;
          success = true;
          e = new_e;

          break;
        } else {
          // Failed to decrease error, increase lambda
          lambda *= lambda_factor;
        }
      }

      result.iterations = i;
      result.H = H;
      result.b = b;
      result.error = e;

      if (!success) {
        break;
      }
    }

    result.num_inliers = std::count_if(factors.begin(), factors.end(), [](const auto& factor) { return factor.inlier(); });

    return result;
  }

  bool verbose;              ///< If true, print debug messages
  int max_iterations;        ///< Max number of optimization iterations
  int max_inner_iterations;  ///< Max number of inner iterations (lambda-trial)
  double init_lambda;        ///< Initial lambda (damping factor)
  double lambda_factor;      ///< Lambda increase factor
};

}  // namespace small_gicp
