// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <small_gicp/util/lie.hpp>
#include <limits>
#include <math/fixed_matrix_ops.h>
#include <small_gicp/factors/linearization.hpp>
#include <small_gicp/ann/traits.hpp>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

/// @brief Point-to-point per-point error factor.
struct ICPFactor {
  struct Setting {};

  ICPFactor(const Setting& setting = Setting()) : target_index(std::numeric_limits<size_t>::max()), source_index(std::numeric_limits<size_t>::max()) {}

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
    const auto& transform = T;
    const auto source_point = traits::point_matrix(source, source_index);
    const auto transformed = transform * source_point;
    const std::array<double, 4> transed_source_pt{{transformed[0][0], transformed[1][0], transformed[2][0], transformed[3][0]}};

    size_t k_index;
    double k_sq_dist;
    if (!traits::nearest_neighbor_search(target_tree, transed_source_pt, &k_index, &k_sq_dist) || rejector(target, source, T, k_index, source_index, k_sq_dist)) {
      return false;
    }

    target_index = k_index;
    const auto target_point = traits::point_matrix(target, target_index);
    const auto residual = target_point + (-1.0) * transformed;

    const auto rotation = rotation_matrix(T);
    const auto jacobian = factor_detail::pointJacobian(
        source_point, rotation);
    factor_detail::linearizeResidual(jacobian, residual, *H, *b, *e);

    return true;
  }

  template <typename TargetPointCloud, typename SourcePointCloud>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T) const {
    if (target_index == std::numeric_limits<size_t>::max()) {
      return 0.0;
    }

    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    const auto& transform = T;
    const auto source_point = traits::point_matrix(source, source_index);
    const auto target_point = traits::point_matrix(target, target_index);
    const auto residual = target_point +
        (-1.0) * (transform * source_point);
    return factor_detail::squaredError(residual);
  }

  bool inlier() const { return target_index != std::numeric_limits<size_t>::max(); }

  size_t target_index;
  size_t source_index;
};
}  // namespace small_gicp
