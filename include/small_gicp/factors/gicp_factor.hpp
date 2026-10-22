// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
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
    using Packed4 = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    using Packed3 = Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
    // Copy Eigen storage only at the point/transform API boundary. All fixed
    // matrix products below use Mantle's existing row-major operations.
    const Packed4 packed_transform = T.matrix();
    const auto transform = Matrix<double, 4, 4>::from(packed_transform.data(), 16);
    const auto source_point = traits::point(source, source_index);
    const auto transformed = transform * Matrix<double, 4, 1>::from(source_point.data(), 4);
    const Eigen::Map<const Eigen::Vector4d> transed_source_pt(transformed.raw_pointer());

    size_t k_index;
    double k_sq_dist;
    if (!traits::nearest_neighbor_search(target_tree, transed_source_pt, &k_index, &k_sq_dist) || rejector(target, source, T, k_index, source_index, k_sq_dist)) {
      return false;
    }

    const Packed4 packed_source_covariance = traits::cov(source, source_index);
    const Packed4 packed_target_covariance = traits::cov(target, k_index);
    const auto source_covariance = Matrix<double, 4, 4>::from(packed_source_covariance.data(), 16);
    const auto target_covariance = Matrix<double, 4, 4>::from(packed_target_covariance.data(), 16);
    Matrix<double, 4, 4> weight{};
    if (!factor_detail::gicpPrecision(target_covariance, source_covariance, transform, weight)) return false;
    mahalanobis = weight;
    target_index = k_index;

    const auto target_point = traits::point(target, target_index);
    const auto residual = Matrix<double, 4, 1>::from(target_point.data(), 4) + (-1.0) * transformed;
    const Packed3 packed_rotation = T.linear();
    const auto rotation = Matrix<double, 3, 3>::from(packed_rotation.data(), 9);
    const auto jacobian = factor_detail::pointJacobian(Matrix<double, 4, 1>::from(source_point.data(), 4), rotation);
    factor_detail::Hessian hessian{};
    factor_detail::Gradient gradient{};
    factor_detail::linearizeResidual(jacobian, residual, weight, hessian, gradient, *e);
    *H = Eigen::Map<const Eigen::Matrix<double, 6, 6, Eigen::RowMajor>>(hessian.raw_pointer());
    *b = Eigen::Map<const Eigen::Matrix<double, 6, 1>>(gradient.raw_pointer());

    return true;
  }

  /// @brief Evaluate error
  /// @param target   Target point cloud
  /// @param source   Source point cloud
  /// @param T        Evaluation point
  /// @return Error
  template <typename TargetPointCloud, typename SourcePointCloud>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Eigen::Isometry3d& T) const {
    if (target_index == std::numeric_limits<size_t>::max()) {
      return 0.0;
    }

    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    using Packed4 = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    const Packed4 packed_transform = T.matrix();
    const auto transform = Matrix<double, 4, 4>::from(packed_transform.data(), 16);
    const auto& weight = mahalanobis;
    const auto source_point = traits::point(source, source_index);
    const auto target_point = traits::point(target, target_index);
    const auto transformed = transform * Matrix<double, 4, 1>::from(source_point.data(), 4);
    const auto residual = Matrix<double, 4, 1>::from(target_point.data(), 4) + (-1.0) * transformed;
    return factor_detail::squaredError(residual, weight);
  }

  /// @brief Returns true if this factor is not rejected as an outlier
  bool inlier() const { return target_index != std::numeric_limits<size_t>::max(); }

  size_t target_index;          ///< Target point index
  size_t source_index;          ///< Source point index
  earth::mantle::collection::Matrix<double, 4, 4> mahalanobis;  ///< Owned fused precision
};
}  // namespace small_gicp
