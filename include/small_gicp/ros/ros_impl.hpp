#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>

#include <small_gicp/points/point_cloud.hpp>

namespace small_gicp {
namespace detail {

template <typename PointFieldT>
struct FieldOffsets {
  std::size_t x = 0;
  std::size_t y = 0;
  std::size_t z = 0;
  bool valid = true;
  bool seen_x = false;
  bool seen_y = false;
  bool seen_z = false;

  bool has_xyz() const { return valid && seen_x && seen_y && seen_z; }

  template <typename FieldsT>
  static FieldOffsets parse(const FieldsT& fields, std::size_t point_step) {
    FieldOffsets offsets;
    for (const auto& field : fields) {
      const bool is_xyz = field.name == "x" || field.name == "y" || field.name == "z";
      if (!is_xyz) continue;
      const bool duplicate = field.name == "x" ? offsets.seen_x : field.name == "y" ? offsets.seen_y : offsets.seen_z;
      const bool layout_ok = field.datatype == PointFieldT::FLOAT32 && field.count == 1 &&
                             static_cast<std::size_t>(field.offset) <= point_step &&
                             sizeof(float) <= point_step - static_cast<std::size_t>(field.offset);
      if (duplicate || !layout_ok) {
        offsets.valid = false;
        continue;
      }
      if (field.name == "x") { offsets.x = field.offset; offsets.seen_x = true; }
      else if (field.name == "y") { offsets.y = field.offset; offsets.seen_y = true; }
      else { offsets.z = field.offset; offsets.seen_z = true; }
    }
    return offsets;
  }
};

inline bool add_overflow(std::size_t a, std::size_t b, std::size_t* out) {
  if (b > std::numeric_limits<std::size_t>::max() - a) return true;
  *out = a + b;
  return false;
}
inline bool mul_overflow(std::size_t a, std::size_t b, std::size_t* out) {
  if (a != 0 && b > std::numeric_limits<std::size_t>::max() / a) return true;
  *out = a * b;
  return false;
}
inline bool host_big_endian() {
  const std::uint16_t value = 0x0102;
  std::uint8_t first = 0;
  std::memcpy(&first, &value, sizeof(first));
  return first == 0x01;
}
inline float read_float(const std::uint8_t* data, bool source_big_endian) {
  std::uint8_t bytes[sizeof(float)];
  std::memcpy(bytes, data, sizeof(bytes));
  if (source_big_endian != host_big_endian()) {
    for (std::size_t i = 0; i < sizeof(bytes) / 2; ++i) std::swap(bytes[i], bytes[sizeof(bytes) - 1 - i]);
  }
  float value = 0.0f;
  std::memcpy(&value, bytes, sizeof(value));
  return value;
}
template <typename PointCloud2T, typename PointFieldT>
PointCloud::Ptr from_impl(const PointCloud2T& msg) {
  auto cloud = std::make_shared<PointCloud>();
  const std::size_t width = msg.width;
  const std::size_t height = msg.height;
  std::size_t num_points = 0;
  if (mul_overflow(width, height, &num_points) || num_points == 0) return cloud;

  const auto offsets = FieldOffsets<PointFieldT>::parse(msg.fields, msg.point_step);
  std::size_t minimum_row_step = 0;
  std::size_t rows_before_last = 0;
  std::size_t last_row_end = 0;
  if (!offsets.has_xyz() || msg.point_step == 0 || mul_overflow(width, msg.point_step, &minimum_row_step) ||
      msg.row_step < minimum_row_step || mul_overflow(height - 1, msg.row_step, &rows_before_last) ||
      add_overflow(rows_before_last, minimum_row_step, &last_row_end) || last_row_end > msg.data.size()) {
    std::cerr << "warning: invalid PointCloud2 XYZ layout" << std::endl;
    return cloud;
  }

  cloud->points.reserve(num_points);
  const std::uint8_t* data_ptr = msg.data.data();
  for (std::size_t row = 0; row < height; ++row) {
    for (std::size_t col = 0; col < width; ++col) {
      const std::uint8_t* pt = data_ptr + row * msg.row_step + col * msg.point_step;
      const float x = read_float(pt + offsets.x, msg.is_bigendian);
      const float y = read_float(pt + offsets.y, msg.is_bigendian);
      const float z = read_float(pt + offsets.z, msg.is_bigendian);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;
      cloud->points.push_back(PointCloud::NativePoint::from({x, y, z, 1.0}));
    }
  }
  cloud->normals.resize(cloud->points.size());
  cloud->covs.resize(cloud->points.size());
  return cloud;
}

template <typename PointCloud2T, typename PointFieldT, typename TimeT>
PointCloud2T to_impl(const PointCloud& cloud, const std::string& frame_id, const TimeT& stamp, bool with_normals) {
  PointCloud2T msg;

  msg.header.frame_id = frame_id;
  msg.header.stamp = stamp;
  msg.height = 1;
  msg.width = static_cast<uint32_t>(cloud.size());
  msg.is_bigendian = false;
  msg.is_dense = true;

  if (cloud.empty()) {
    msg.point_step = 0;
    msg.row_step = 0;
    return msg;
  }

  // Build fields
  PointFieldT field;
  field.count = 1;
  field.datatype = PointFieldT::FLOAT32;

  uint32_t offset = 0;

  field.name = "x";
  field.offset = offset;
  msg.fields.push_back(field);
  offset += sizeof(float);

  field.name = "y";
  field.offset = offset;
  msg.fields.push_back(field);
  offset += sizeof(float);

  field.name = "z";
  field.offset = offset;
  msg.fields.push_back(field);
  offset += sizeof(float);

  if (with_normals) {
    field.name = "normal_x";
    field.offset = offset;
    msg.fields.push_back(field);
    offset += sizeof(float);

    field.name = "normal_y";
    field.offset = offset;
    msg.fields.push_back(field);
    offset += sizeof(float);

    field.name = "normal_z";
    field.offset = offset;
    msg.fields.push_back(field);
    offset += sizeof(float);
  }

  msg.point_step = offset;
  msg.row_step = msg.point_step * msg.width;

  // Fill data
  msg.data.resize(msg.row_step);
  uint8_t* data_ptr = msg.data.data();

  for (size_t i = 0; i < cloud.size(); i++) {
    uint8_t* pt = data_ptr + i * msg.point_step;

    const auto& point = cloud.point(i);
    *reinterpret_cast<float*>(pt + 0) = static_cast<float>(point[0][0]);
    *reinterpret_cast<float*>(pt + 4) = static_cast<float>(point[1][0]);
    *reinterpret_cast<float*>(pt + 8) = static_cast<float>(point[2][0]);

    if (with_normals) {
      const auto& normal = cloud.normal(i);
      *reinterpret_cast<float*>(pt + 12) = static_cast<float>(normal[0][0]);
      *reinterpret_cast<float*>(pt + 16) = static_cast<float>(normal[1][0]);
      *reinterpret_cast<float*>(pt + 20) = static_cast<float>(normal[2][0]);
    }
  }

  return msg;
}

}  // namespace detail

}  // namespace small_gicp
