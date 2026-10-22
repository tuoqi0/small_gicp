// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <cstddef>
#include <iterator>
#include <collection/fixed_matrix.h>
#include <math/symmetric_eigen3.h>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

/// @brief Parameters to control the projection axis search.
struct ProjectionSetting {
  int max_scan_count = 128;  ///< Maximum number of points to use for the axis search.
};

/// @brief Conventional axis-aligned projection (i.e., selecting any of XYZ axes with the largest variance).
/// @note  Up to max_scan_count samples are used to estimate the variance.
struct AxisAlignedProjection {
public:
  /// @brief Project the point to the selected axis.
  /// @param pt  Point to project
  /// @return    Projected value
  template <typename Point>
  double operator()(const Point& pt) const { return pt[axis]; }

  /// @brief Find the axis with the largest variance.
  /// @param points     Point cloud
  /// @param first      First point index iterator
  /// @param last       Last point index iterator
  /// @param setting    Search setting
  /// @return           Projection with the largest variance axis
  template <typename PointCloud, typename IndexConstIterator>
  static AxisAlignedProjection find_axis(const PointCloud& points, IndexConstIterator first, IndexConstIterator last, const ProjectionSetting& setting) {
    const size_t N = std::distance(first, last);
    earth::mantle::collection::Matrix<double, 4, 1> sum_pt{};
    earth::mantle::collection::Matrix<double, 3, 1> sum_sq{};

    const size_t max_scan_count = static_cast<size_t>(setting.max_scan_count);
    const size_t step = N < max_scan_count ? 1 : N / max_scan_count;
    const size_t num_steps = N / step;
    for (size_t i = 0; i < num_steps; i++) {
      const auto itr = first + step * i;
      const auto pt = traits::point(points, *itr);
      for (size_t axis = 0; axis < 3; ++axis) {
        sum_pt[axis][0] += pt[axis];
        sum_sq[axis][0] += pt[axis] * pt[axis];
      }
      sum_pt[3][0] += pt[3];
    }

    std::array<double, 3> var{};
    for (size_t axis = 0; axis < 3; ++axis)
      var[axis] = sum_sq[axis][0] - (sum_pt[axis][0] / sum_pt[3][0]) * sum_pt[axis][0];

    return AxisAlignedProjection{var[0] > var[1] ? (var[0] > var[2] ? 0 : 2) : (var[1] > var[2] ? 1 : 2)};
  }

public:
  int axis;  ///< Axis index (0: X, 1: Y, 2: Z)
};

/// @brief Normal projection (i.e., selecting the 3D direction with the largest variance of the points).
/// @note  Up to max_scan_count samples are used to estimate the variance along the axis.
struct NormalProjection {
public:
  /// @brief Project the point to the normal direction.
  /// @param pt   Point to project
  /// @return     Projected value
  template <typename Point>
  double operator()(const Point& pt) const {
    return normal[0] * pt[0] + normal[1] * pt[1] + normal[2] * pt[2];
  }

  /// @brief  Find the direction with the largest variance.
  /// @param points   Point cloud
  /// @param first    First point index iterator
  /// @param last     Last point index iterator
  /// @param setting  Search setting
  /// @return         Projection with the largest variance direction
  template <typename PointCloud, typename IndexConstIterator>
  static NormalProjection find_axis(const PointCloud& points, IndexConstIterator first, IndexConstIterator last, const ProjectionSetting& setting) {
    const size_t N = std::distance(first, last);
    NormalProjection p{{1, 0, 0}};
    earth::mantle::collection::Matrix<double, 4, 1> sum_pt{};
    earth::mantle::collection::Matrix<double, 3, 3> sum_sq{};

    const size_t max_scan_count = static_cast<size_t>(setting.max_scan_count);
    const size_t step = N < max_scan_count ? 1 : N / max_scan_count;
    const size_t num_steps = N / step;
    for (size_t i = 0; i < num_steps; i++) {
      const auto itr = first + step * i;
      const auto pt = traits::point(points, *itr);
      for (size_t row = 0; row < 3; ++row) {
        sum_pt[row][0] += pt[row];
        for (size_t col = 0; col < 3; ++col)
          sum_sq[row][col] += pt[row] * pt[col];
      }
      sum_pt[3][0] += pt[3];
    }

    std::array<double, 9> covariance{};
    for (std::size_t row = 0; row < 3; ++row) {
      const double mean = sum_pt[row][0] / sum_pt[3][0];
      for (std::size_t col = 0; col < 3; ++col)
        covariance[row * 3 + col] = (sum_sq[row][col] - mean * sum_pt[col][0]) / sum_pt[3][0];
    }
    earth::mantle::math::SymmetricEigen3 spectrum;
    // A coordinate axis is still a valid KD-tree split if covariance is invalid.
    if (!earth::mantle::math::symmetricEigen3(covariance, spectrum)) return p;
    for (std::size_t row = 0; row < 3; ++row) p.normal[row] = spectrum.vectors[row * 3 + 2];
    return p;
  }

public:
  std::array<double, 3> normal;  ///< Projection direction
};

}  // namespace small_gicp
