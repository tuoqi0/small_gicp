// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>
#include <collection/fixed_matrix.h>
#include <math/symmetric_eigen3.h>
#include <small_gicp/ann/kdtree.hpp>

namespace small_gicp {

/// @brief Computes point normals from eigenvectors and sets them to the point cloud.
/// @note  If a sufficient number of neighbor points is not found, a zero vector is set to the point.
template <typename PointCloud>
struct NormalSetter {
  /// @brief Handle invalid case (too few points).
  static void set_invalid(PointCloud& cloud, size_t i) { traits::set_normal(cloud, i, std::array<double, 4>{}); }

  /// @brief Compute and set the normal to the point cloud.
  static void set(PointCloud& cloud, size_t i, const earth::mantle::math::SymmetricEigen3& eig) {
    const double x = eig.vectors[0], y = eig.vectors[3], z = eig.vectors[6];
    const double norm = std::sqrt(x * x + y * y + z * z);
    std::array<double, 4> normal{x / norm, y / norm, z / norm, 0.0};
    const auto point = traits::point(cloud, i);
    const double facing = (traits::detail::coordinate(point, 0) * normal[0] + traits::detail::coordinate(point, 2) * normal[2]) +
                          (traits::detail::coordinate(point, 1) * normal[1] + traits::detail::coordinate(point, 3) * normal[3]);
    if (facing > 0) for (double& value : normal) value = -value;
    traits::set_normal(cloud, i, normal);
  }
};

/// @brief Computes point covariances from eigenvectors and sets them to the point cloud.
/// @note  If a sufficient number of neighbor points is not found, an identity matrix is set to the point.
template <typename PointCloud>
struct CovarianceSetter {
  /// @brief Handle invalid case (too few points).
  static void set_invalid(PointCloud& cloud, size_t i) {
    earth::mantle::collection::Matrix<double, 4, 4> cov{};
    for (size_t axis = 0; axis < 3; ++axis) cov[axis][axis] = 1.0;
    traits::set_cov(cloud, i, cov);
  }

  /// @brief Compute and set the covariance to the point cloud.
  static void set(PointCloud& cloud, size_t i, const earth::mantle::math::SymmetricEigen3& eig) {
    const std::array<double, 3> values{{1e-3, 1.0, 1.0}};
    earth::mantle::collection::Matrix<double, 4, 4> cov{};
    for (size_t row = 0; row < 3; ++row) {
      for (size_t col = 0; col < 3; ++col) {
        for (size_t axis = 0; axis < 3; ++axis)
          cov[row][col] += eig.vectors[row * 3 + axis] * values[axis] * eig.vectors[col * 3 + axis];
      }
    }
    traits::set_cov(cloud, i, cov);
  }
};

/// @brief Computes point normals and covariances from eigenvectors and sets them to the point cloud.
/// @note  If a sufficient number of neighbor points is not found, a zero vector and an identity matrix are set to the point.
template <typename PointCloud>
struct NormalCovarianceSetter {
  /// @brief Handle invalid case (too few points).
  static void set_invalid(PointCloud& cloud, size_t i) {
    NormalSetter<PointCloud>::set_invalid(cloud, i);
    CovarianceSetter<PointCloud>::set_invalid(cloud, i);
  }

  /// @brief Compute and set the normal and covariance to the point cloud.
  static void set(PointCloud& cloud, size_t i, const earth::mantle::math::SymmetricEigen3& eig) {
    NormalSetter<PointCloud>::set(cloud, i, eig);
    CovarianceSetter<PointCloud>::set(cloud, i, eig);
  }
};

template <typename Setter, typename PointCloud, typename Tree>
void estimate_local_features(PointCloud& cloud, Tree& kdtree, int num_neighbors, size_t point_index) {
  std::vector<size_t> k_indices(num_neighbors);
  std::vector<double> k_sq_dists(num_neighbors);
  const size_t n = kdtree.knn_search(traits::point(cloud, point_index), num_neighbors, k_indices.data(), k_sq_dists.data());

  if (n < 5) {
    // Insufficient number of neighbors
    Setter::set_invalid(cloud, point_index);
    return;
  }

  earth::mantle::collection::Matrix<double, 3, 1> sum_points{};
  earth::mantle::collection::Matrix<double, 3, 3> sum_cross{};
  for (size_t i = 0; i < n; i++) {
    const auto& pt = traits::point(cloud, k_indices[i]);
    for (size_t row = 0; row < 3; ++row) {
      sum_points[row][0] += traits::detail::coordinate(pt, row);
      for (size_t col = 0; col <= row; ++col)
        sum_cross[row][col] += traits::detail::coordinate(pt, row) * traits::detail::coordinate(pt, col);
    }
  }

  // Match the original population covariance and its lower-triangle input.
  // The homogeneous fourth coordinate does not enter the 3D eigensystem.
  std::array<double, 9> covariance{};
  const double count = static_cast<double>(n);
  for (size_t row = 0; row < 3; ++row) {
    const double mean = sum_points[row][0] / count;
    for (size_t col = 0; col <= row; ++col) {
      const double value = (sum_cross[row][col] - mean * sum_points[col][0]) / count;
      covariance[row * 3 + col] = value;
      covariance[col * 3 + row] = value;
    }
  }
  earth::mantle::math::SymmetricEigen3 eig;
  if (!earth::mantle::math::symmetricEigen3(covariance, eig)) {
    Setter::set_invalid(cloud, point_index);
    return;
  }
  Setter::set(cloud, point_index, eig);
}

template <typename Setter, typename PointCloud, typename KdTree>
void estimate_local_features(PointCloud& cloud, const KdTree& kdtree, int num_neighbors) {
  traits::resize(cloud, traits::size(cloud));
  for (size_t i = 0; i < traits::size(cloud); i++) {
    estimate_local_features<Setter>(cloud, kdtree, num_neighbors, i);
  }
}

template <typename Setter, typename PointCloud>
void estimate_local_features(PointCloud& cloud, int num_neighbors) {
  traits::resize(cloud, traits::size(cloud));

  UnsafeKdTree<PointCloud> kdtree(cloud);
  for (size_t i = 0; i < traits::size(cloud); i++) {
    estimate_local_features<Setter>(cloud, kdtree, num_neighbors, i);
  }
}

/// @brief Estimate point normals
/// @note  If a sufficient number of neighbor points for normal/covariance estimation (5 points) is not found,
///        an invalid normal/covariance is set to the point (normal=zero vector, covariance=identity matrix).
/// @param cloud          [in/out] Point cloud
/// @param num_neighbors  Number of neighbors used for attribute estimation
template <typename PointCloud>
void estimate_normals(PointCloud& cloud, int num_neighbors = 20) {
  estimate_local_features<NormalSetter<PointCloud>>(cloud, num_neighbors);
}

/// @brief Estimate point normals
/// @note  If a sufficient number of neighbor points for normal/covariance estimation (5 points) is not found,
///        an invalid normal/covariance is set to the point (normal=zero vector, covariance=identity matrix).
/// @param cloud          [in/out] Point cloud
/// @param kdtree         Nearest neighbor search
/// @param num_neighbors  Number of neighbors used for attribute estimation
template <typename PointCloud, typename KdTree>
void estimate_normals(PointCloud& cloud, KdTree& kdtree, int num_neighbors = 20) {
  estimate_local_features<NormalSetter<PointCloud>>(cloud, kdtree, num_neighbors);
}

/// @brief Estimate point covariances
/// @note  If a sufficient number of neighbor points for normal/covariance estimation (5 points) is not found,
///        an invalid normal/covariance is set to the point (normal=zero vector, covariance=identity matrix).
/// @param cloud          [in/out] Point cloud
/// @param num_neighbors  Number of neighbors used for attribute estimation
template <typename PointCloud>
void estimate_covariances(PointCloud& cloud, int num_neighbors = 20) {
  estimate_local_features<CovarianceSetter<PointCloud>>(cloud, num_neighbors);
}

/// @brief Estimate point covariances
/// @note  If a sufficient number of neighbor points for normal/covariance estimation (5 points) is not found,
///        an invalid normal/covariance is set to the point (normal=zero vector, covariance=identity matrix).
/// @param cloud          [in/out] Point cloud
/// @param kdtree         Nearest neighbor search
/// @param num_neighbors  Number of neighbors used for attribute estimation
template <typename PointCloud, typename KdTree>
void estimate_covariances(PointCloud& cloud, KdTree& kdtree, int num_neighbors = 20) {
  estimate_local_features<CovarianceSetter<PointCloud>>(cloud, kdtree, num_neighbors);
}

/// @brief Estimate point normals and covariances
/// @note  If a sufficient number of neighbor points for normal/covariance estimation (5 points) is not found,
///        an invalid normal/covariance is set to the point (normal=zero vector, covariance=identity matrix).
/// @param cloud          [in/out] Point cloud
/// @param num_neighbors  Number of neighbors used for attribute estimation
template <typename PointCloud>
void estimate_normals_covariances(PointCloud& cloud, int num_neighbors = 20) {
  estimate_local_features<NormalCovarianceSetter<PointCloud>>(cloud, num_neighbors);
}

/// @brief Estimate point normals and covariances
/// @note  If a sufficient number of neighbor points for normal/covariance estimation (5 points) is not found,
///        an invalid normal/covariance is set to the point (normal=zero vector, covariance=identity matrix).
/// @param cloud          [in/out] Point cloud
/// @param kdtree         Nearest neighbor search
/// @param num_neighbors  Number of neighbors used for attribute estimation
template <typename PointCloud, typename KdTree>
void estimate_normals_covariances(PointCloud& cloud, KdTree& kdtree, int num_neighbors = 20) {
  estimate_local_features<NormalCovarianceSetter<PointCloud>>(cloud, kdtree, num_neighbors);
}

}  // namespace small_gicp
