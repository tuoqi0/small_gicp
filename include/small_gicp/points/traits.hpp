// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <collection/fixed_matrix.h>

namespace small_gicp {

namespace traits {

template <typename T>
struct Traits;

/// @brief  Get the number of points.
template <typename T>
size_t size(const T& points) {
  return Traits<T>::size(points);
}

/// @brief Check if the point cloud has points.
template <typename T>
bool has_points(const T& points) {
  return Traits<T>::has_points(points);
}

/// @brief Check if the point cloud has normals.
template <typename T>
bool has_normals(const T& points) {
  return Traits<T>::has_normals(points);
}

/// @brief Check if the point cloud has covariances.
template <typename T>
bool has_covs(const T& points) {
  return Traits<T>::has_covs(points);
}

/// @brief Get i-th point. 4D vector is used to take advantage of SIMD intrinsics. The last element must be filled by one (x, y, z, 1).
template <typename T>
auto point(const T& points, size_t i) {
  return Traits<T>::point(points, i);
}

/// @brief Get i-th normal. 4D vector is used to take advantage of SIMD intrinsics. The last element must be filled by zero (nx, ny, nz, 0).
template <typename T>
auto normal(const T& points, size_t i) {
  return Traits<T>::normal(points, i);
}

/// @brief Get i-th covariance. Only the top-left 3x3 matrix is filled, and the bottom row and the right col must be filled by zero.
template <typename T>
auto cov(const T& points, size_t i) {
  return Traits<T>::cov(points, i);
}

namespace detail {
template <typename Point>
auto coordinate(const Point& point, size_t axis) -> decltype(point[axis]) { return point[axis]; }
template <typename T, size_t Dimension>
const T& coordinate(const earth::mantle::collection::Matrix<T, Dimension, 1>& point, size_t axis) {
  return point[axis][0];
}
template <typename Covariance>
auto coefficient(const Covariance& covariance, size_t row, size_t col) -> decltype(covariance(row, col)) {
  return covariance(row, col);
}
template <typename T, size_t Rows, size_t Cols>
const T& coefficient(const earth::mantle::collection::Matrix<T, Rows, Cols>& covariance, size_t row, size_t col) {
  return covariance[row][col];
}
}  // namespace detail

// Preserve direct homogeneous division used by the existing voxel centroids.
inline earth::mantle::collection::Matrix<double, 4, 1> normalized_point(
    const earth::mantle::collection::Matrix<double, 4, 1>& sum) {
  earth::mantle::collection::Matrix<double, 4, 1> result{};
  for (size_t axis = 0; axis < 4; ++axis) result[axis][0] = sum[axis][0] / sum[3][0];
  return result;
}

// Explicit point-cloud boundary. Native matrices pass by value, as the
// original generic getters did; external layouts are read coordinate-wise.
template <typename T>
earth::mantle::collection::Matrix<double, 4, 1> point_matrix(const T& points, size_t i) {
  const auto value = point(points, i);
  return earth::mantle::collection::Matrix<double, 4, 1>::from({
      detail::coordinate(value, 0), detail::coordinate(value, 1),
      detail::coordinate(value, 2), detail::coordinate(value, 3)});
}
template <typename T>
earth::mantle::collection::Matrix<double, 4, 1> normal_matrix(const T& points, size_t i) {
  const auto value = normal(points, i);
  return earth::mantle::collection::Matrix<double, 4, 1>::from({
      detail::coordinate(value, 0), detail::coordinate(value, 1),
      detail::coordinate(value, 2), detail::coordinate(value, 3)});
}
template <typename T>
earth::mantle::collection::Matrix<double, 4, 4> covariance_matrix(const T& points, size_t i) {
  const auto value = cov(points, i);
  earth::mantle::collection::Matrix<double, 4, 4> result{};
  for (size_t row = 0; row < 4; ++row)
    for (size_t col = 0; col < 4; ++col) result[row][col] = detail::coefficient(value, row, col);
  return result;
}

/// @brief Resize the point cloud (this function should resize all attributes)
template <typename T>
void resize(T& points, size_t n) {
  Traits<T>::resize(points, n);
}

/// @brief Set i-th point. (x, y, z, 1)
template <typename T, typename Point>
void set_point(T& points, size_t i, const Point& pt) {
  Traits<T>::set_point(points, i, pt);
}

/// @brief Set i-th normal. (nx, nz, nz, 0)
template <typename T, typename Normal>
void set_normal(T& points, size_t i, const Normal& pt) {
  Traits<T>::set_normal(points, i, pt);
}

/// @brief Set i-th covariance. Only the top-left 3x3 matrix should be filled.
template <typename T, typename Covariance>
void set_cov(T& points, size_t i, const Covariance& cov) {
  Traits<T>::set_cov(points, i, cov);
}

}  // namespace traits
}  // namespace small_gicp
