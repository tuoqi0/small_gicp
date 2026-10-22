// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <small_gicp/util/lie.hpp>
#include <math/fixed_matrix_ops.h>
#include <math/positive_definite_solve.h>
#include <small_gicp/factors/linearization.hpp>
#include <small_gicp/ann/traits.hpp>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

/// @brief GICP (distribution-to-distribution) per-point error factor.
struct GICPFactor {
  struct Setting {};

  /// @brief Constructor
  GICPFactor(const Setting& setting = Setting())
  : target_index(std::numeric_limits<size_t>::max()),
    source_index(std::numeric_limits<size_t>::max()),
    mahalanobis{} {}

  /// @brief Linearize the factor
  /// @param target       Target point cloud
  /// @param source       Source point cloud
  /// @param target_tree  Nearest neighbor search for the target point cloud
  /// @param T            Linearization point
  /// @param source_index Source point index
  /// @param rejector     Correspondence rejector
  /// @param H            Linearized information matrix
  /// @param b            Linearized information vector
  /// @param e            Error at the linearization point
  /// @return             True if the point is inlier
  template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree, typename CorrespondenceRejector>
  bool linearize(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const Transform& T,
    size_t source_index,
    const CorrespondenceRejector& rejector,
    factor_detail::Hessian* H,
    factor_detail::Gradient* b,
    double* e) {
    //
    this->source_index = source_index;
    this->target_index = std::numeric_limits<size_t>::max();

    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    // Point-cloud layout conversion is isolated in traits; transforms are native.
    const auto& transform = T;
    const auto source_point = traits::point_matrix(source, source_index);
    const auto transformed = transform * source_point;
    const std::array<double, 4> transed_source_pt{{transformed[0][0], transformed[1][0], transformed[2][0], transformed[3][0]}};

    size_t k_index;
    double k_sq_dist;
    if (!traits::nearest_neighbor_search(target_tree, transed_source_pt, &k_index, &k_sq_dist) || rejector(target, source, T, k_index, source_index, k_sq_dist)) {
      return false;
    }

    const auto source_covariance = traits::covariance_matrix(source, source_index);
    const auto target_covariance = traits::covariance_matrix(target, k_index);
    Matrix<double, 4, 4> weight{};
    if (!factor_detail::gicpPrecision(target_covariance, source_covariance, transform, weight)) return false;
    mahalanobis = weight;
    target_index = k_index;

    const auto target_point = traits::point_matrix(target, target_index);
    const auto residual = target_point + (-1.0) * transformed;
    const auto rotation = rotation_matrix(T);
    const auto jacobian = factor_detail::pointJacobian(source_point, rotation);
    factor_detail::linearizeResidual(jacobian, residual, weight, *H, *b, *e);

    return true;
  }

  /// @brief Evaluate error
  /// @param target   Target point cloud
  /// @param source   Source point cloud
  /// @param T        Evaluation point
  /// @return Error
  template <typename TargetPointCloud, typename SourcePointCloud>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T) const {
    if (target_index == std::numeric_limits<size_t>::max()) {
      return 0.0;
    }

    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    const auto& transform = T;
    const auto& weight = mahalanobis;
    const auto source_point = traits::point_matrix(source, source_index);
    const auto target_point = traits::point_matrix(target, target_index);
    const auto transformed = transform * source_point;
    const auto residual = target_point + (-1.0) * transformed;
    return factor_detail::squaredError(residual, weight);
  }

  /// @brief Returns true if this factor is not rejected as an outlier
  bool inlier() const { return target_index != std::numeric_limits<size_t>::max(); }

  size_t target_index;          ///< Target point index
  size_t source_index;          ///< Source point index
  earth::mantle::collection::Matrix<double, 4, 4> mahalanobis;  ///< Owned fused precision
};
}  // namespace small_gicp
