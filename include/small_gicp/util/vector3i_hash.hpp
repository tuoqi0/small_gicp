// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <cstddef>

namespace small_gicp {

using VoxelCoord = std::array<int, 3>;

/**
 * @brief Spatial hashing function.
 *        Teschner et al., "Optimized Spatial Hashing for Collision Detection of Deformable Objects", VMV2003.
 */
struct XORVector3iHash {
public:
  template <typename Coordinate>
  size_t operator()(const Coordinate& x) const {
    const size_t p1 = 73856093;
    const size_t p2 = 19349669;  // 19349663 was not a prime number
    const size_t p3 = 83492791;
    return static_cast<size_t>((x[0] * p1) ^ (x[1] * p2) ^ (x[2] * p3));
  }

  template <typename Coordinate>
  static size_t hash(const Coordinate& x) { return XORVector3iHash()(x); }
  template <typename Coordinate>
  static bool equal(const Coordinate& x1, const Coordinate& x2) { return x1 == x2; }
};

}  // namespace small_gicp
