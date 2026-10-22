// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <collection/fixed_matrix.h>
#include <small_gicp/util/lie.hpp>

namespace small_gicp {

// Explicit interoperability boundary. Views borrow storage; never map a temporary.
template <typename T, size_t Rows, size_t Cols>
using EigenMatrix = Eigen::Matrix<T, Rows, Cols, Cols == 1 ? Eigen::ColMajor : Eigen::RowMajor>;

template <typename T, size_t Rows, size_t Cols>
Eigen::Map<EigenMatrix<T, Rows, Cols>> eigen_view(earth::mantle::collection::Matrix<T, Rows, Cols>& value) {
  return Eigen::Map<EigenMatrix<T, Rows, Cols>>(value.raw_pointer());
}
template <typename T, size_t Rows, size_t Cols>
Eigen::Map<const EigenMatrix<T, Rows, Cols>> eigen_view(const earth::mantle::collection::Matrix<T, Rows, Cols>& value) {
  return Eigen::Map<const EigenMatrix<T, Rows, Cols>>(value.raw_pointer());
}
template <typename T, size_t Rows, size_t Cols>
void eigen_view(const earth::mantle::collection::Matrix<T, Rows, Cols>&&) = delete;

template <typename Derived>
auto native_matrix(const Eigen::MatrixBase<Derived>& value) {
  static_assert(Derived::RowsAtCompileTime > 0 && Derived::ColsAtCompileTime > 0,
                "fixed-size matrix required");
  earth::mantle::collection::Matrix<typename Derived::Scalar,
      Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> result{};
  eigen_view(result) = value;
  return result;
}

inline Transform native_transform(const Eigen::Isometry3d& value) {
  return native_matrix(value.matrix());
}
inline Eigen::Isometry3d eigen_transform(const Transform& value) {
  return Eigen::Isometry3d(eigen_view(value));
}
inline Eigen::Matrix3d skew(const Eigen::Vector3d& value) {
  const auto native = native_matrix(value);
  const auto result = earth::mantle::math::hat(native);
  return eigen_view(result);
}
inline Eigen::Quaterniond so3_exp(const Eigen::Vector3d& value) {
  const auto q = earth::mantle::math::so3Exp({{value[0], value[1], value[2]}});
  return Eigen::Quaterniond(q[3], q[0], q[1], q[2]);
}
inline Eigen::Isometry3d se3_exp(const Eigen::Matrix<double, 6, 1>& value) {
  return eigen_transform(se3_exp(native_matrix(value)));
}

}  // namespace small_gicp
