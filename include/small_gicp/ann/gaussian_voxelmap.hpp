// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math/fixed_matrix_ops.h>

#include <small_gicp/ann/traits.hpp>
#include <small_gicp/points/traits.hpp>
#include <small_gicp/ann/incremental_voxelmap.hpp>

namespace small_gicp {

/// @brief Gaussian voxel that computes and stores voxel mean and covariance.
struct GaussianVoxel {
public:
  struct Setting {};

  /// @brief Constructor.
  GaussianVoxel() : finalized(false), num_points(0), mean{}, cov{} {}

  /// @brief  Number of points in the voxel (Always 1 for GaussianVoxel).
  size_t size() const { return 1; }

  /// @brief  Add a point to the voxel.
  /// @param  setting        Setting
  /// @param  transformed_pt Transformed point mean
  /// @param  points         Point cloud
  /// @param  i              Index of the point
  /// @param  T              Transformation matrix
  template <typename PointCloud>
  void add(const Setting& setting, const Eigen::Vector4d& transformed_pt, const PointCloud& points, size_t i, const Eigen::Isometry3d& T) {
    using earth::mantle::math::operator*;
    using earth::mantle::math::operator+;
    using earth::mantle::math::transpose;
    using PackedCovariance = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    const PackedCovariance packed_transform = T.matrix();
    const PackedCovariance packed_covariance = traits::cov(points, i);
    const auto transform = NativeCovariance::from(packed_transform.data(), 16);
    const auto covariance = NativeCovariance::from(packed_covariance.data(), 16);
    if (finalized) {
      this->finalized = false;
      this->mean = static_cast<double>(num_points) * this->mean;
      this->cov = static_cast<double>(num_points) * this->cov;
    }

    num_points++;
    this->mean = this->mean + NativePoint::from(transformed_pt.data(), 4);
    this->cov = this->cov + transform * covariance * transpose(transform);
  }

  /// @brief Finalize the voxel mean and covariance.
  void finalize() {
    if (finalized) {
      return;
    }

    // Preserve direct division (not reciprocal multiplication) and the original
    // sum/average lifecycle when more points arrive after finalization.
    const double count = static_cast<double>(num_points);
    for (size_t row = 0; row < 4; ++row) {
      mean[row][0] /= count;
      for (size_t col = 0; col < 4; ++col) cov[row][col] /= count;
    }
    finalized = true;
  }

public:
  bool finalized;        ///< If true, mean and cov are finalized, otherwise they represent the sum of input points
  size_t num_points;     ///< Number of input points
  using NativePoint = earth::mantle::collection::Matrix<double, 4, 1>;
  using NativeCovariance = earth::mantle::collection::Matrix<double, 4, 4>;
  NativePoint mean;      ///< Mean (homogeneous coordinate included)
  NativeCovariance cov;  ///< Row-major covariance
};

namespace traits {

template <>
struct Traits<GaussianVoxel> {
  static size_t size(const GaussianVoxel& voxel) { return 1; }
  static bool has_points(const GaussianVoxel& voxel) { return true; }
  static bool has_covs(const GaussianVoxel& voxel) { return true; }

  static Eigen::Vector4d point(const GaussianVoxel& voxel, size_t i) {
    return Eigen::Map<const Eigen::Vector4d>(voxel.mean.raw_pointer());
  }
  static Eigen::Matrix4d cov(const GaussianVoxel& voxel, size_t i) {
    using PackedCovariance = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    return Eigen::Map<const PackedCovariance>(voxel.cov.raw_pointer());
  }

  static size_t nearest_neighbor_search(const GaussianVoxel& voxel, const Eigen::Vector4d& pt, size_t* k_index, double* k_sq_dist) {
    *k_index = 0;
    *k_sq_dist = (point(voxel, 0) - pt).squaredNorm();
    return 1;
  }

  static size_t knn_search(const GaussianVoxel& voxel, const Eigen::Vector4d& pt, size_t k, size_t* k_index, double* k_sq_dist) {
    return nearest_neighbor_search(voxel, pt, k_index, k_sq_dist);
  }

  template <typename Result>
  static void knn_search(const GaussianVoxel& voxel, const Eigen::Vector4d& pt, Result& result) {
    result.push(0, (point(voxel, 0) - pt).squaredNorm());
  }
};

}  // namespace traits

using GaussianVoxelMap = IncrementalVoxelMap<GaussianVoxel>;

}  // namespace small_gicp
