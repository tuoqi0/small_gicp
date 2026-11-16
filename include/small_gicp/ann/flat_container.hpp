#pragma once

#include <algorithm>
#include <type_traits>
#include <vector>
#include <small_gicp/util/lie.hpp>
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
  void add(const Setting& setting, const NativePoint& transformed_pt, const PointCloud& points, size_t i, const Transform& T) {
    if (this->points.size() >= setting.max_num_points_in_cell) return;
    if (std::any_of(this->points.begin(), this->points.end(), [&](const auto& pt) {
        return ann_detail::squared_distance(pt, transformed_pt) < setting.min_sq_dist_in_cell;
      })) {
      return;
    }

    using earth::mantle::math::operator*;
    using earth::mantle::math::transpose;
    this->points.push_back(transformed_pt);
    if constexpr (HasNormals || HasCovs) {
      const auto& transform = T;
      if constexpr (HasNormals) {
        const auto normal = traits::normal_matrix(points, i);
        this->normals.push_back(transform * normal);
      }
      if constexpr (HasCovs) {
        const auto covariance = traits::covariance_matrix(points, i);
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
  template <typename Query>
  size_t nearest_neighbor_search(const Query& pt, size_t* k_index, double* k_sq_dist) const {
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
  template <typename Query>
  size_t knn_search(const Query& pt, int k, size_t* k_indices, double* k_sq_dists) const {
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
  template <typename Query, typename Result>
  void knn_search(const Query& pt, Result& result) const {
    if (points.empty()) {
      return;
    }

    for (size_t i = 0; i < points.size(); i++) {
      const double sq_dist = ann_detail::squared_distance(points[i], pt);
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

  static auto point(const FlatContainer<HasNormals, HasCovs>& container, size_t i) { return container.points[i]; }
  static auto normal(const FlatContainer<HasNormals, HasCovs>& container, size_t i) { return container.normals[i]; }
  static auto cov(const FlatContainer<HasNormals, HasCovs>& container, size_t i) { return container.covs[i]; }

  template <typename Query>
  static size_t nearest_neighbor_search(const FlatContainer<HasNormals, HasCovs>& container, const Query& pt, size_t* k_index, double* k_sq_dist) {
    return container.nearest_neighbor_search(pt, k_index, k_sq_dist);
  }

  template <typename Query>
  static size_t knn_search(const FlatContainer<HasNormals, HasCovs>& container, const Query& pt, size_t k, size_t* k_index, double* k_sq_dist) {
    return container.knn_search(pt, k, k_index, k_sq_dist);
  }

  template <typename Query, typename Result>
  static void knn_search(const FlatContainer<HasNormals, HasCovs>& container, const Query& pt, Result& result) {
    container.knn_search(pt, result);
  }
};

}  // namespace traits

}  // namespace small_gicp
