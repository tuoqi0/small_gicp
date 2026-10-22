// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <small_gicp/points/traits.hpp>
#include <small_gicp/util/vector3i_hash.hpp>

namespace small_gicp {

namespace detail {
inline int floor_coordinate(double value) {
  const int truncated = static_cast<int>(value);
  return truncated - (value < static_cast<double>(truncated));
}
}

/// @brief Fast floor (https://stackoverflow.com/questions/824118/why-is-floor-so-slow).
/// @param pt  Double vector
/// @return    Floored int vector
template <typename Point>
std::array<int, 4> fast_floor(const Point& pt) {
  std::array<int, 4> coord{};
  for (size_t axis = 0; axis < coord.size(); ++axis) {
    coord[axis] = detail::floor_coordinate(traits::detail::coordinate(pt, axis));
  }
  return coord;
}

template <typename Point>
VoxelCoord voxel_coordinates(const Point& pt, double scale) {
  VoxelCoord coord{};
  for (size_t axis = 0; axis < coord.size(); ++axis) {
    coord[axis] = detail::floor_coordinate(traits::detail::coordinate(pt, axis) * scale);
  }
  return coord;
}

}  // namespace small_gicp
