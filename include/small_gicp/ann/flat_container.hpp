// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <type_traits>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math/fixed_matrix_ops.h>
#include <small_gicp/ann/traits.hpp>
#include <small_gicp/points/traits.hpp>
#include <small_gicp/ann/knn_result.hpp>

namespace small_gicp {

/// @brief Point container with a flat vector.
/// @note  IncrementalVoxelMap combined with FlastContainer is mostly the same as linear iVox.
///        Bai et al., "Faster-LIO: Lightweight Tightly Coupled Lidar-Inertial Odometry Using Parallel Sparse Incremental Voxels", IEEE RA-L, 2022
/// @note  This container stores only up to max_num_points_in_cell points and avoids insertings points that are too close to existing points (min_sq_dist_in_cell).
/// @tparam HasNormals  If true, store normals.
/// @tparam HasCovs     If true, store covariances.
template <bool HasNormals = false, bool HasCovs = false>
struct FlatContainer {
public:
  using NativePoint = earth::mantle::collection::Matrix<double, 4, 1>;
  using NativeCovariance = earth::mantle::collection::Matrix<double, 4, 4>;
  /// @brief FlatContainer setting.
  struct Setting {
    double min_sq_dist_in_cell = 0.1 * 0.1;  ///< Minimum squared distance between points in a cell.
    size_t max_num_points_in_cell = 10;      ///< Maximum number of points in a cell.
  };

  /// @brief Constructor.
  FlatContainer() { points.reserve(5); }

  /// @brief Number of points.
  size_t size() const { return points.size(); }

  /// @brief Add a point to the container.
  ///        Ignore the input when the cell is full or an existing point is too close.
  /// @param setting         Setting
  /// @param transformed_pt  Transformed point (== T * points[i])
  /// @param points          Point cloud
  /// @param i               Index of the point
  /// @param T               Transformation matrix
  template <typename PointCloud>
  void add(const Setting& setting, const Eigen::Vector4d& transformed_pt, const PointCloud& points, size_t i, const Eigen::Isometry3d& T) {
    if (this->points.size() >= setting.max_num_points_in_cell) return;
    if (std::any_of(this->points.begin(), this->points.end(), [&](const auto& pt) {
        return (Eigen::Map<const Eigen::Vector4d>(pt.raw_pointer()) - transformed_pt).squaredNorm() < setting.min_sq_dist_in_cell;
      })) {
      return;
    }

    using earth::mantle::math::operator*;
    using earth::mantle::math::transpose;
    this->points.push_back(NativePoint::from(transformed_pt.data(), 4));
    if constexpr (HasNormals || HasCovs) {
      using PackedCovariance = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
      const PackedCovariance packed_transform = T.matrix();
      const auto transform = NativeCovariance::from(packed_transform.data(), 16);
      if constexpr (HasNormals) {
        const auto normal = traits::normal(points, i);
        this->normals.push_back(transform * NativePoint::from(normal.data(), 4));
      }
      if constexpr (HasCovs) {
        const PackedCovariance packed_covariance = traits::cov(points, i);
        const auto covariance = NativeCovariance::from(packed_covariance.data(), 16);
        this->covs.push_back(transform * covariance * transpose(transform));
      }
    }
  }

  /// @brief Finalize the container (Nothing to do for FlatContainer).
  void finalize() {}

  /// @brief Find the nearest neighbor.
  /// @param pt           Query point
  /// @param k_index      Index of the nearest neighbor
  /// @param k_sq_dist    Squared distance to the nearest neighbor
  /// @return             Number of found points (0 or 1)
  size_t nearest_neighbor_search(const Eigen::Vector4d& pt, size_t* k_index, double* k_sq_dist) const {
    if (points.empty()) {
      return 0;
    }

    KnnResult<1> result(k_index, k_sq_dist);
    knn_search(pt, result);
    return result.num_found();
  }

  /// @brief Find k nearest neighbors.
  /// @param pt           Query point
  /// @param k            Number of neighbors
  /// @param k_indices    Indices of nearest neighbors
  /// @param k_sq_dists   Squared distances to nearest neighbors (sorted in ascending order)
  /// @return             Number of found points
  size_t knn_search(const Eigen::Vector4d& pt, int k, size_t* k_indices, double* k_sq_dists) const {
    if (points.empty()) {
      return 0;
    }

    KnnResult<-1> result(k_indices, k_sq_dists, k);
    knn_search(pt, result);
    return result.num_found();
  }

  /// @brief Find k nearest neighbors and push them into the given result.
  /// @param pt           Query point
  /// @param result       Result collector that receives the candidates
  template <typename Result>
  void knn_search(const Eigen::Vector4d& pt, Result& result) const {
    if (points.empty()) {
      return;
    }

    for (size_t i = 0; i < points.size(); i++) {
      const double sq_dist = (Eigen::Map<const Eigen::Vector4d>(points[i].raw_pointer()) - pt).squaredNorm();
      result.push(i, sq_dist);
    }
  }

public:
  struct Empty {};

  std::vector<NativePoint> points;                                          ///< Points
  std::conditional_t<HasNormals, std::vector<NativePoint>, Empty> normals;  ///< Normals (Empty if HasNormals is false)
  std::conditional_t<HasCovs, std::vector<NativeCovariance>, Empty> covs;    ///< Covariances (Empty if HasCovs is false)
};

/// @brief FlatContainer that stores only points.
using FlatContainerPoints = FlatContainer<false, false>;
/// @brief FlatContainer with normals.
using FlatContainerNormal = FlatContainer<true, false>;
/// @brief FlatContainer with covariances.
using FlatContainerCov = FlatContainer<false, true>;
/// @brief FlatContainer with normals and covariances.
using FlatContainerNormalCov = FlatContainer<true, true>;

namespace traits {

template <bool HasNormals, bool HasCovs>
struct Traits<FlatContainer<HasNormals, HasCovs>> {
  static size_t size(const FlatContainer<HasNormals, HasCovs>& container) { return container.size(); }
  static bool has_points(const FlatContainer<HasNormals, HasCovs>& container) { return container.size(); }
  static bool has_normals(const FlatContainer<HasNormals, HasCovs>& container) { return HasNormals && container.size(); }
  static bool has_covs(const FlatContainer<HasNormals, HasCovs>& container) { return HasCovs && container.size(); }

  static Eigen::Vector4d point(const FlatContainer<HasNormals, HasCovs>& container, size_t i) {
    return Eigen::Map<const Eigen::Vector4d>(container.points[i].raw_pointer());
  }
  static Eigen::Vector4d normal(const FlatContainer<HasNormals, HasCovs>& container, size_t i) {
    return Eigen::Map<const Eigen::Vector4d>(container.normals[i].raw_pointer());
  }
  static Eigen::Matrix4d cov(const FlatContainer<HasNormals, HasCovs>& container, size_t i) {
    using PackedCovariance = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    return Eigen::Map<const PackedCovariance>(container.covs[i].raw_pointer());
  }

  static size_t nearest_neighbor_search(const FlatContainer<HasNormals, HasCovs>& container, const Eigen::Vector4d& pt, size_t* k_index, double* k_sq_dist) {
    return container.nearest_neighbor_search(pt, k_index, k_sq_dist);
  }

  static size_t knn_search(const FlatContainer<HasNormals, HasCovs>& container, const Eigen::Vector4d& pt, size_t k, size_t* k_index, double* k_sq_dist) {
    return container.knn_search(pt, k, k_index, k_sq_dist);
  }

  template <typename Result>
  static void knn_search(const FlatContainer<HasNormals, HasCovs>& container, const Eigen::Vector4d& pt, Result& result) {
    container.knn_search(pt, result);
  }
};

}  // namespace traits

}  // namespace small_gicp
