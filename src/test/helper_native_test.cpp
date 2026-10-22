#include <robot_core/test/doctest.h>
#include <small_gicp/registration/registration_helper.hpp>
#include <array>
#include <vector>

namespace {
std::vector<std::array<double, 3>> makeCorner() {
  std::vector<std::array<double, 3>> points;
  for (int row = -5; row <= 5; ++row) {
    for (int col = -5; col <= 5; ++col) {
      points.push_back({{4, row * .25, col * .25}});
      points.push_back({{row * .25, 3, col * .25}});
      points.push_back({{row * .25, col * .25, -1}});
    }
  }
  return points;
}

template <typename Point>
void checkRawInput() {
  const auto coordinates = makeCorner();
  std::vector<Point> points(coordinates.size());
  for (size_t i = 0; i < points.size(); ++i)
    for (size_t axis = 0; axis < 3; ++axis) points[i][axis] = coordinates[i][axis];
  const auto prepared = small_gicp::preprocess_points(points, .1, 10, 1);
  REQUIRE(prepared.first->size() == points.size());
  CHECK(prepared.second != nullptr);
  CHECK(prepared.first->point(0)[3][0] == 1);
}
}

TEST_CASE("public helper preprocesses standard coordinate arrays without Eigen") {
  checkRawInput<std::array<float, 3>>();
  checkRawInput<std::array<float, 4>>();
  checkRawInput<std::array<double, 3>>();
  checkRawInput<std::array<double, 4>>();
  const small_gicp::PointCloud native(makeCorner());
  const auto prepared = small_gicp::preprocess_points(native.points, .1, 10, 1);
  CHECK(prepared.first->size() == native.size());
  const auto empty = small_gicp::preprocess_points(std::vector<std::array<double, 3>>{}, .1, 10, 1);
  CHECK(empty.first->empty());
  CHECK(small_gicp::create_gaussian_voxelmap(*empty.first, .5)->size() == 0);
}

TEST_CASE("public helper preserves native initial transform for every registration mode") {
  const auto target = makeCorner();
  auto source = target;
  auto expected = earth::mantle::math::identityMatrix<double, 4>();
  expected[0][3] = .08;
  expected[1][3] = -.06;
  expected[2][3] = .04;
  for (auto& point : source)
    for (size_t axis = 0; axis < 3; ++axis) point[axis] -= expected[axis][3];
  for (const int threads : {1, 2}) {
    for (const auto type : {small_gicp::RegistrationSetting::ICP,
        small_gicp::RegistrationSetting::PLANE_ICP,
        small_gicp::RegistrationSetting::GICP,
        small_gicp::RegistrationSetting::VGICP}) {
      small_gicp::RegistrationSetting setting;
      setting.type = type;
      setting.num_threads = threads;
      setting.downsampling_resolution = .1;
      setting.voxel_resolution = .2;
      const auto result = small_gicp::align(target, source, expected, setting);
      REQUIRE(result.converged);
      CHECK(result.num_inliers == source.size());
      for (size_t row = 0; row < 4; ++row)
        for (size_t col = 0; col < 4; ++col)
          CHECK(result.T_target_source[row][col] == doctest::Approx(expected[row][col]).epsilon(1e-8));
    }
  }
  const auto result = small_gicp::align(target, target);
  CHECK(result.converged);
  CHECK(result.num_inliers == target.size());
}
