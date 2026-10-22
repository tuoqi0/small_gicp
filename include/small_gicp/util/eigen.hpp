// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <Eigen/Core>
#include <collection/fixed_matrix.h>

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

}  // namespace small_gicp
