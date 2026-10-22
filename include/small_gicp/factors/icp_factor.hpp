// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
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
    const Eigen::Isometry3d& T,
    size_t source_index,
    const CorrespondenceRejector& rejector,
    Eigen::Matrix<double, 6, 6>* H,
    Eigen::Matrix<double, 6, 1>* b,
    double* e) {
    //
    this->source_index = source_index;
    this->target_index = std::numeric_limits<size_t>::max();

    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    using PackedTransform = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    using PackedRotation = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
    const PackedTransform packed_transform = T.matrix();
    const auto transform = Matrix<double, 4, 4>::from(packed_transform.data(), packed_transform.size());
    const auto source_point = traits::point(source, source_index);
    const auto transformed = transform * Matrix<double, 4, 1>::from(source_point.data(), source_point.size());
    const Eigen::Map<const Eigen::Vector4d> transed_source_pt(transformed.raw_pointer());

    size_t k_index;
    double k_sq_dist;
    if (!traits::nearest_neighbor_search(target_tree, transed_source_pt, &k_index, &k_sq_dist) || rejector(target, source, T, k_index, source_index, k_sq_dist)) {
      return false;
    }

    target_index = k_index;
    const auto target_point = traits::point(target, target_index);
    const auto residual = Matrix<double, 4, 1>::from(target_point.data(), target_point.size()) + (-1.0) * transformed;

    const PackedRotation packed_rotation = T.linear();
    const auto rotation = Matrix<double, 3, 3>::from(packed_rotation.data(), packed_rotation.size());
    const auto jacobian = factor_detail::pointJacobian(
        Matrix<double, 4, 1>::from(source_point.data(), source_point.size()), rotation);
    factor_detail::Hessian hessian{};
    factor_detail::Gradient gradient{};
    factor_detail::linearizeResidual(jacobian, residual, hessian, gradient, *e);
    *H = Eigen::Map<const Eigen::Matrix<double, 6, 6, Eigen::RowMajor>>(hessian.raw_pointer());
    *b = Eigen::Map<const Eigen::Matrix<double, 6, 1>>(gradient.raw_pointer());

    return true;
  }

  template <typename TargetPointCloud, typename SourcePointCloud>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Eigen::Isometry3d& T) const {
    if (target_index == std::numeric_limits<size_t>::max()) {
      return 0.0;
    }

    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    using PackedTransform = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    const PackedTransform packed_transform = T.matrix();
    const auto transform = Matrix<double, 4, 4>::from(packed_transform.data(), packed_transform.size());
    const auto source_point = traits::point(source, source_index);
    const auto target_point = traits::point(target, target_index);
    const auto residual = Matrix<double, 4, 1>::from(target_point.data(), target_point.size()) +
        (-1.0) * (transform * Matrix<double, 4, 1>::from(source_point.data(), source_point.size()));
    return factor_detail::squaredError(residual);
  }

  bool inlier() const { return target_index != std::numeric_limits<size_t>::max(); }

  size_t target_index;
  size_t source_index;
};
}  // namespace small_gicp
