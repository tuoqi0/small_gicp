#pragma once

#include <small_gicp/util/lie.hpp>
#include <limits>
#include <math/fixed_matrix_ops.h>
#include <small_gicp/factors/linearization.hpp>
#include <small_gicp/ann/traits.hpp>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

/// @brief Point-to-plane per-point error factor.
struct PointToPlaneICPFactor {
  struct Setting {};

  PointToPlaneICPFactor(const Setting& setting = Setting()) : target_index(std::numeric_limits<size_t>::max()), source_index(std::numeric_limits<size_t>::max()) {}

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
    const auto target_normal = traits::normal_matrix(target, target_index);

    const auto target_point = traits::point_matrix(target, target_index);
    const auto residual = target_point + (-1.0) * transformed;
    Matrix<double, 4, 1> weighted_residual{};
    for (size_t row = 0; row < 4; ++row)
      weighted_residual[row][0] = target_normal[row][0] * residual[row][0];

    // Preserve diag(normal) * R * skew(point), not a scalar plane projection.
    auto weighted_rotation = rotation_matrix(T);
    for (size_t row = 0; row < 3; ++row)
      for (size_t col = 0; col < 3; ++col) weighted_rotation[row][col] *= target_normal[row][0];
    const auto jacobian = factor_detail::pointJacobian(
        source_point, weighted_rotation);
    factor_detail::linearizeResidual(jacobian, weighted_residual, *H, *b, *e);

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
    const auto normal = traits::normal_matrix(target, target_index);
    const auto residual = target_point +
        (-1.0) * (transform * source_point);
    Matrix<double, 4, 1> weighted_residual{};
    for (size_t row = 0; row < 4; ++row)
      weighted_residual[row][0] = normal[row][0] * residual[row][0];
    return factor_detail::squaredError(weighted_residual);
  }

  bool inlier() const { return target_index != std::numeric_limits<size_t>::max(); }

  size_t target_index;
  size_t source_index;
};
}  // namespace small_gicp
