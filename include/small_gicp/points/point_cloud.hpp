// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <memory>
#include <array>
#include <vector>
#include <Eigen/Core>
#include <collection/fixed_matrix.h>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

/**
 * @brief Point cloud
 */
struct PointCloud {
public:
  using Ptr = std::shared_ptr<PointCloud>;
  using ConstPtr = std::shared_ptr<const PointCloud>;
  using NativePoint = earth::mantle::collection::Matrix<double, 4, 1>;
  using NativeCovariance = earth::mantle::collection::Matrix<double, 4, 4>;
  using PointView = Eigen::Map<Eigen::Vector4d>;
  using ConstPointView = Eigen::Map<const Eigen::Vector4d>;
  using PackedCovariance = Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
  using CovarianceView = Eigen::Map<PackedCovariance>;
  using ConstCovarianceView = Eigen::Map<const PackedCovariance>;

  /// @brief Constructor
  PointCloud() {}
  PointCloud(const PointCloud&) = default;
  PointCloud& operator=(const PointCloud&) = default;
  PointCloud(PointCloud&&) noexcept = default;
  PointCloud& operator=(PointCloud&&) noexcept = default;

  /// @brief Constructor
  /// @param points  Points to initialize the point cloud
  template <typename T, int D, typename Allocator>
  explicit PointCloud(const std::vector<Eigen::Matrix<T, D, 1>, Allocator>& points) {
    this->resize(points.size());
    for (size_t i = 0; i < points.size(); i++) {
      this->point(i) << points[i].template cast<double>().template head<3>(), 1.0;
    }
  }

  /// @brief Destructor
  ~PointCloud() {}

  /// @brief Number of points.
  size_t size() const { return points.size(); }

  /// @brief Check if the point cloud is empty.
  bool empty() const { return points.empty(); }

  /// @brief Resize point/normal/cov buffers.
  /// @param n  Number of points
  void resize(size_t n) {
    points.resize(n);
    normals.resize(n);
    covs.resize(n);
  }

  /// @brief Get i-th point.
  PointView point(size_t i) { return PointView(points[i].raw_pointer()); }

  /// @brief Get i-th normal.
  PointView normal(size_t i) { return PointView(normals[i].raw_pointer()); }

  /// @brief Get i-th covariance.
  CovarianceView cov(size_t i) { return CovarianceView(covs[i].raw_pointer()); }

  /// @brief Get i-th point (const).
  ConstPointView point(size_t i) const {
    return ConstPointView(points[i].raw_pointer());
  }

  /// @brief Get i-th normal (const).
  ConstPointView normal(size_t i) const {
    return ConstPointView(normals[i].raw_pointer());
  }

  /// @brief Get i-th covariance (const).
  ConstCovarianceView cov(size_t i) const {
    return ConstCovarianceView(covs[i].raw_pointer());
  }

public:
  // Native buffers own all data. Getter views borrow storage and, like the
  // previous references, must not survive resize/destruction of their buffer.
  // Generic traits below return values to retain their previous copy semantics.
  std::vector<NativePoint> points;      ///< Point coordinates (x, y, z, 1)
  std::vector<NativePoint> normals;     ///< Point normals (nx, ny, nz, 0)
  std::vector<NativeCovariance> covs;    ///< Row-major covariance with zero padding
};

namespace traits {

template <>
struct Traits<PointCloud> {
  using Points = PointCloud;

  static size_t size(const Points& points) { return points.size(); }

  static bool has_points(const Points& points) { return !points.points.empty(); }
  static bool has_normals(const Points& points) { return !points.normals.empty(); }
  static bool has_covs(const Points& points) { return !points.covs.empty(); }

  static Eigen::Vector4d point(const Points& points, size_t i) { return points.point(i); }
  static Eigen::Vector4d normal(const Points& points, size_t i) { return points.normal(i); }
  static Eigen::Matrix4d cov(const Points& points, size_t i) { return points.cov(i); }

  static void resize(Points& points, size_t n) { points.resize(n); }
  static void set_point(Points& points, size_t i, const Eigen::Vector4d& pt) { points.point(i) = pt; }
  static void set_normal(Points& points, size_t i, const Eigen::Vector4d& n) { points.normal(i) = n; }
  static void set_cov(Points& points, size_t i, const Eigen::Matrix4d& cov) { points.cov(i) = cov; }
  static void set_normal(Points& points, size_t i, const std::array<double, 4>& normal) {
    points.normal(i) = Eigen::Map<const Eigen::Vector4d>(normal.data());
  }
  static void set_cov(Points& points, size_t i, const Points::NativeCovariance& cov) {
    points.covs[i] = cov;
  }
};

}  // namespace traits

}  // namespace small_gicp
