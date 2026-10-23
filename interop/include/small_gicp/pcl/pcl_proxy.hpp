// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <pcl/point_cloud.h>
#include <collection/fixed_matrix.h>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

/// @brief Proxy class to access PCL point cloud with external covariance matrices.
template <typename PointT>
struct PointCloudProxy {
  PointCloudProxy(const pcl::PointCloud<PointT>& cloud, std::vector<Eigen::Matrix4d>& covs) : cloud(cloud), covs(covs) {}

  const pcl::PointCloud<PointT>& cloud; ///< Point cloud
  std::vector<Eigen::Matrix4d>& covs;   ///< Covariance matrices
};

namespace traits {
template <typename PointT>
struct Traits<PointCloudProxy<PointT>> {
  using Points = PointCloudProxy<PointT>;

  static size_t size(const Points& points) { return points.cloud.size(); }

  static bool has_points(const Points& points) { return !points.cloud.points.empty(); }
  static bool has_covs(const Points& points) { return !points.covs.empty(); }

  static const Eigen::Vector4d point(const Points& points, size_t i) { return points.cloud.at(i).getVector4fMap().template cast<double>(); }
  static const Eigen::Matrix4d& cov(const Points& points, size_t i) { return points.covs[i]; }

  static void resize(Points& points, size_t n) { points.covs.resize(n); }
  static void set_cov(Points& points, size_t i, const Eigen::Matrix4d& cov) { points.covs[i] = cov; }
  static void set_cov(Points& points, size_t i, const earth::mantle::collection::Matrix<double, 4, 4>& cov) {
    using RowMajorCovariance = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
    points.covs[i] = Eigen::Map<const RowMajorCovariance>(cov.raw_pointer());
  }
};

}  // namespace traits

}  // namespace small_gicp
