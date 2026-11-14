#include <robot_core/test/doctest.h>
#include <small_gicp/ann/projection.hpp>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

struct NativeProjectionPoints {
  std::vector<std::array<double, 4>> points;
};

namespace small_gicp { namespace traits {
template <>
struct Traits<NativeProjectionPoints> {
  static const std::array<double, 4>& point(const NativeProjectionPoints& cloud, size_t index) {
    return cloud.points[index];
  }
};
}}

TEST_CASE("projection accepts native points without Eigen or PCL headers") {
  for (const double scale : {1e-6, 1.0, 1e6}) {
    NativeProjectionPoints cloud;
    for (int i = -30; i <= 30; ++i) {
      const double t = scale * i;
      cloud.points.push_back({{t, 2 * t, 3 * t, 1}});
    }
    std::vector<size_t> indices(cloud.points.size());
    std::iota(indices.begin(), indices.end(), 0);
    small_gicp::ProjectionSetting setting;
    setting.max_scan_count = 9;
    const auto aligned = small_gicp::AxisAlignedProjection::find_axis(
      cloud, indices.begin(), indices.end(), setting);
    const auto normal = small_gicp::NormalProjection::find_axis(
      cloud, indices.begin(), indices.end(), setting);
    CHECK(aligned.axis == 2);
    const double alignment = normal.normal[0] + 2 * normal.normal[1] + 3 * normal.normal[2];
    CHECK(std::abs(alignment) == doctest::Approx(std::sqrt(14)).epsilon(1e-12));
    const std::array<double, 4> query{{2, 4, 6, 99}};
    CHECK(aligned(query) == 6);
    CHECK(normal(query) == doctest::Approx(2 * alignment).epsilon(1e-12));
  }
}

TEST_CASE("axis projection retains the original variance tie order") {
  NativeProjectionPoints cloud{{{{-1, -1, 0, 1}}, {{1, 1, 0, 1}}}};
  const std::array<size_t, 2> indices{{0, 1}};
  const auto equal_xy = small_gicp::AxisAlignedProjection::find_axis(
    cloud, indices.begin(), indices.end(), small_gicp::ProjectionSetting{});
  CHECK(equal_xy.axis == 1);
  cloud.points[0][2] = -1;
  cloud.points[1][2] = 1;
  const auto equal_xyz = small_gicp::AxisAlignedProjection::find_axis(
    cloud, indices.begin(), indices.end(), small_gicp::ProjectionSetting{});
  CHECK(equal_xyz.axis == 2);
}

TEST_CASE("native normal projection retains degenerate and nonfinite covariance handling") {
  NativeProjectionPoints cloud{{{{1, 2, 3, 1}}}};
  const std::array<size_t, 1> indices{{0}};
  const auto single_point = small_gicp::NormalProjection::find_axis(
    cloud, indices.begin(), indices.end(), small_gicp::ProjectionSetting{});
  CHECK(single_point.normal == std::array<double, 3>{{0, 0, 1}});
  cloud.points[0][0] = std::numeric_limits<double>::quiet_NaN();
  const auto nonfinite = small_gicp::NormalProjection::find_axis(
    cloud, indices.begin(), indices.end(), small_gicp::ProjectionSetting{});
  CHECK(nonfinite.normal == std::array<double, 3>{{1, 0, 0}});
}
