// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <atomic>
#include <memory>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <numeric>
#include <math/fixed_matrix_ops.h>
#include <iostream>

#include <tbb/tbb.h>
#include <small_gicp/points/traits.hpp>
#include <small_gicp/util/fast_floor.hpp>
#include <small_gicp/util/vector3i_hash.hpp>

namespace small_gicp {

/// @brief Voxel grid downsampling with TBB backend.
/// @note  This function has minor run-by-run non-deterministic behavior due to parallel data collection that results
///        in a deviation of the number of points in the downsampling results (up to 10% increase from the single-thread version).
/// @note  Discretized voxel coords must be in 21bit range [-1048576, 1048575].
///        For example, if the downsampling resolution is 0.01 m, point coordinates must be in [-10485.76, 10485.75] m.
///        Points outside the valid range will be ignored.
/// @param points     Input points
/// @param leaf_size  Downsampling resolution
/// @return           Downsampled points
template <typename InputPointCloud, typename OutputPointCloud = InputPointCloud>
std::shared_ptr<OutputPointCloud> voxelgrid_sampling_tbb(const InputPointCloud& points, double leaf_size) {
  using earth::mantle::math::operator*;
  using earth::mantle::math::operator+;
  if (traits::size(points) == 0) {
    return std::make_shared<OutputPointCloud>();
  }

  const double inv_leaf_size = 1.0 / leaf_size;

  constexpr std::uint64_t invalid_coord = std::numeric_limits<std::uint64_t>::max();
  constexpr int coord_bit_size = 21;                       // Bits to represent each voxel coordinate (pack 21x3 = 63bits in 64bit int)
  constexpr size_t coord_bit_mask = (1 << 21) - 1;         // Bit mask
  constexpr int coord_offset = 1 << (coord_bit_size - 1);  // Coordinate offset to make values positive

  std::vector<std::pair<std::uint64_t, size_t>> coord_pt(traits::size(points));
  tbb::parallel_for(tbb::blocked_range<size_t>(0, traits::size(points), 64), [&](const tbb::blocked_range<size_t>& range) {
    for (size_t i = range.begin(); i != range.end(); i++) {
      auto coord = fast_floor(traits::point_matrix(points, i) * inv_leaf_size);
      for (int& value : coord) value += coord_offset;
      if (std::any_of(coord.begin(), coord.end(), [=](int value) { return value < 0 || value > static_cast<int>(coord_bit_mask); })) {
        std::cerr << "warning: voxel coord is out of range!!" << std::endl;
        coord_pt[i] = {invalid_coord, i};
        continue;
      }

      // Compute voxel coord bits (0|1bit, z|21bit, y|21bit, x|21bit)
      const std::uint64_t bits =                                                           //
        (static_cast<std::uint64_t>(coord[0] & coord_bit_mask) << (coord_bit_size * 0)) |  //
        (static_cast<std::uint64_t>(coord[1] & coord_bit_mask) << (coord_bit_size * 1)) |  //
        (static_cast<std::uint64_t>(coord[2] & coord_bit_mask) << (coord_bit_size * 2));
      coord_pt[i] = {bits, i};
    }
  });

  // Sort by voxel coords
  tbb::parallel_sort(coord_pt, [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

  auto downsampled = std::make_shared<OutputPointCloud>();
  traits::resize(*downsampled, traits::size(points));

  // Take block-wise sum
  const int block_size = 2048;
  std::atomic_uint64_t num_points = 0;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, traits::size(points), block_size), [&](const tbb::blocked_range<size_t>& range) {
    std::vector<earth::mantle::collection::Matrix<double, 4, 1>> sub_points;
    sub_points.reserve(block_size);

    earth::mantle::collection::Matrix<double, 4, 1> sum_pt = traits::point_matrix(points, coord_pt[range.begin()].second);
    for (size_t i = range.begin() + 1; i != range.end(); i++) {
      if (coord_pt[i].first == invalid_coord) {
        continue;
      }

      if (coord_pt[i - 1].first != coord_pt[i].first) {
        sub_points.emplace_back(traits::normalized_point(sum_pt));
        sum_pt = {};
      }
      sum_pt = sum_pt + traits::point_matrix(points, coord_pt[i].second);
    }
    sub_points.emplace_back(traits::normalized_point(sum_pt));

    const size_t point_index_begin = num_points.fetch_add(sub_points.size());
    for (size_t i = 0; i < sub_points.size(); i++) {
      traits::set_point(*downsampled, point_index_begin + i, sub_points[i]);
    }
  });

  traits::resize(*downsampled, num_points);

  return downsampled;
}

}  // namespace small_gicp
