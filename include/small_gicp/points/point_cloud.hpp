// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <memory>
#include <array>
#include <vector>
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

  /// @brief Constructor
  PointCloud() {}
  PointCloud(const PointCloud&) = default;
  PointCloud& operator=(const PointCloud&) = default;
  PointCloud(PointCloud&&) noexcept = default;
  PointCloud& operator=(PointCloud&&) noexcept = default;

  /// @brief Constructor
  /// @param points  Points to initialize the point cloud
  template <typename Point, typename Allocator>
  explicit PointCloud(const std::vector<Point, Allocator>& points) {
    this->resize(points.size());
    for (size_t i = 0; i < points.size(); i++) {
      this->points[i] = NativePoint::from({static_cast<double>(traits::detail::coordinate(points[i], 0)),
          static_cast<double>(traits::detail::coordinate(points[i], 1)),
          static_cast<double>(traits::detail::coordinate(points[i], 2)), 1.0});
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
  NativePoint& point(size_t i) { return points[i]; }

  /// @brief Get i-th normal.
  NativePoint& normal(size_t i) { return normals[i]; }

  /// @brief Get i-th covariance.
  NativeCovariance& cov(size_t i) { return covs[i]; }

  /// @brief Get i-th point (const).
  const NativePoint& point(size_t i) const {
    return points[i];
  }

  /// @brief Get i-th normal (const).
  const NativePoint& normal(size_t i) const {
    return normals[i];
  }

  /// @brief Get i-th covariance (const).
  const NativeCovariance& cov(size_t i) const {
    return covs[i];
  }

public:
  // Native buffers own all data. Getter references must not survive resize or
  // destruction of their buffer; no third-party matrix view is stored here.
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

  static Points::NativePoint point(const Points& points, size_t i) { return points.point(i); }
  static Points::NativePoint normal(const Points& points, size_t i) { return points.normal(i); }
  static Points::NativeCovariance cov(const Points& points, size_t i) { return points.cov(i); }

  static void resize(Points& points, size_t n) { points.resize(n); }
  template <typename Point>
  static void set_point(Points& points, size_t i, const Point& pt) {
    points.point(i) = Points::NativePoint::from({detail::coordinate(pt, 0), detail::coordinate(pt, 1),
        detail::coordinate(pt, 2), detail::coordinate(pt, 3)});
  }
  template <typename Point>
  static void set_normal(Points& points, size_t i, const Point& normal) {
    points.normal(i) = Points::NativePoint::from({detail::coordinate(normal, 0), detail::coordinate(normal, 1),
        detail::coordinate(normal, 2), detail::coordinate(normal, 3)});
  }
  template <typename Covariance>
  static void set_cov(Points& points, size_t i, const Covariance& cov) {
    Points::NativeCovariance copy{};
    for (size_t row = 0; row < 4; ++row)
      for (size_t col = 0; col < 4; ++col) copy[row][col] = detail::coefficient(cov, row, col);
    points.covs[i] = copy;
  }
};

}  // namespace traits

}  // namespace small_gicp
