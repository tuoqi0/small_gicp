#include <robot_core/test/doctest.h>
#include <small_gicp/ann/gaussian_voxelmap.hpp>
#include <small_gicp/ann/sequential_voxelmap_accessor.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#include <array>
#include <cmath>
#include <type_traits>
#include <vector>

namespace {
using Voxel = small_gicp::GaussianVoxel;
using Cloud = small_gicp::PointCloud;
using Point = Cloud::NativePoint;
static_assert(std::is_same<decltype(Voxel::mean), Point>::value, "voxel mean must own native storage");
static_assert(std::is_same<decltype(Voxel::cov), Cloud::NativeCovariance>::value, "voxel covariance must own native storage");
struct AcceptCorrespondence {
  template <typename... Args> bool operator()(const Args&...) const { return false; }
};
void checkPoint(const Point& point, std::array<double, 4> expected) {
  for (size_t axis = 0; axis < 4; ++axis) CHECK(point[axis][0] == doctest::Approx(expected[axis]));
}
}

TEST_CASE("native Gaussian voxel preserves analytic transformed accumulation and finalize lifecycle") {
  Cloud cloud; cloud.resize(1);
  cloud.point(0) = Point::from({1, 2, 3, 1});
  cloud.cov(0)[0][0] = 1; cloud.cov(0)[1][1] = 2; cloud.cov(0)[2][2] = 4;
  auto transform = earth::mantle::math::identityMatrix<double, 4>();
  transform[0][0] = transform[1][1] = 0;
  transform[0][1] = -1; transform[1][0] = 1;
  transform[0][3] = 10; transform[1][3] = -2; transform[2][3] = 1;
  Voxel voxel;
  voxel.add(Voxel::Setting{}, Point::from({8, -1, 4, 1}), cloud, 0, transform);
  voxel.finalize(); voxel.finalize();
  REQUIRE(voxel.num_points == 1); CHECK(voxel.finalized);
  checkPoint(voxel.mean, {{8, -1, 4, 1}});
  for (size_t row = 0; row < 4; ++row)
    for (size_t col = 0; col < 4; ++col)
      CHECK(voxel.cov[row][col] == (row == col ? std::array<double, 4>{{2, 1, 4, 0}}[row] : 0));
  cloud.point(0) = Point::from({3, 4, 5, 1});
  voxel.add(Voxel::Setting{}, Point::from({6, 1, 6, 1}), cloud, 0, transform);
  CHECK_FALSE(voxel.finalized); REQUIRE(voxel.num_points == 2);
  voxel.finalize(); checkPoint(voxel.mean, {{7, 0, 5, 1}});
  CHECK(voxel.cov[0][0] == 2); CHECK(voxel.cov[1][1] == 1); CHECK(voxel.cov[2][2] == 4);
}

TEST_CASE("native Gaussian voxel preserves covariance layout and value snapshots") {
  Cloud cloud; cloud.resize(1); cloud.point(0) = Point::from({1, 2, 3, 1});
  for (int row = 0; row < 4; ++row)
    for (int col = 0; col < 4; ++col) cloud.cov(0)[row][col] = 10 * row + col;
  const auto identity = earth::mantle::math::identityMatrix<double, 4>();
  Voxel voxel; voxel.add(Voxel::Setting{}, cloud.point(0), cloud, 0, identity); voxel.finalize();
  CHECK(voxel.cov[0][1] == 1); CHECK(voxel.cov[1][0] == 10);
  const auto point_snapshot = small_gicp::traits::point(voxel, 0);
  const auto covariance_snapshot = small_gicp::traits::cov(voxel, 0);
  Voxel copy = voxel; copy.mean[0][0] = 9; CHECK(voxel.mean[0][0] == 1);
  cloud.point(0) = Point::from({3, 4, 5, 1});
  voxel.add(Voxel::Setting{}, cloud.point(0), cloud, 0, identity); voxel.finalize();
  checkPoint(point_snapshot, {{1, 2, 3, 1}}); checkPoint(voxel.mean, {{2, 3, 4, 1}});
  for (size_t row = 0; row < 4; ++row)
    for (size_t col = 0; col < 4; ++col) CHECK(covariance_snapshot[row][col] == cloud.cov(0)[row][col]);
}

TEST_CASE("native Gaussian map supports repeated insertion search and GICP linearization") {
  Cloud cloud(std::vector<std::array<double, 3>>{{{.1, .2, .3}}});
  cloud.cov(0) = earth::mantle::math::identityMatrix<double, 4>(); cloud.cov(0)[3][3] = 0;
  small_gicp::GaussianVoxelMap map(1.0); map.insert(cloud);
  cloud.point(0)[0][0] = .3; map.insert(cloud);
  REQUIRE(map.size() == 1);
  const small_gicp::SequentialVoxelMapAccessor<small_gicp::GaussianVoxelMap> accessor(map);
  REQUIRE(small_gicp::traits::size(accessor) == 1);
  checkPoint(small_gicp::traits::point(accessor, 0), {{.2, .2, .3, 1}});
  size_t index = 0; double distance = 0;
  REQUIRE(map.nearest_neighbor_search(cloud.point(0), &index, &distance) == 1);
  CHECK(distance == doctest::Approx(.01));
  small_gicp::GICPFactor factor;
  small_gicp::factor_detail::Hessian hessian{};
  small_gicp::factor_detail::Gradient gradient{};
  double error = 0; const auto identity = earth::mantle::math::identityMatrix<double, 4>();
  REQUIRE(factor.linearize(map, cloud, map, identity, 0, AcceptCorrespondence{}, &hessian, &gradient, &error));
  CHECK(error == doctest::Approx(.0025)); CHECK(factor.error(map, cloud, identity) == doctest::Approx(error));
}

TEST_CASE("native Gaussian map registration recovers a rigid pose") {
  using earth::mantle::math::operator*;
  std::vector<std::array<double, 4>> points, scan;
  for (int row = -6; row <= 6; ++row) for (int col = -6; col <= 6; ++col) {
    points.push_back({{4, row * .25, col * .25, 1}});
    points.push_back({{row * .25, 3, col * .25, 1}});
    points.push_back({{row * .25, col * .25, -1, 1}});
  }
  const double angle = .015 / std::sqrt(6.);
  auto expected = small_gicp::se3_exp(small_gicp::factor_detail::Gradient::from({angle, 2 * angle, -angle, 0, 0, 0}));
  expected[0][3] = .08; expected[1][3] = -.06; expected[2][3] = .04;
  const auto inverse_rotation = earth::mantle::math::transpose(small_gicp::rotation_matrix(expected));
  for (const auto& point : points) {
    const auto shifted = earth::mantle::collection::Matrix<double, 3, 1>::from({point[0] - .08, point[1] + .06, point[2] - .04});
    const auto transformed = inverse_rotation * shifted;
    scan.push_back({{transformed[0][0], transformed[1][0], transformed[2][0], 1}});
  }
  Cloud target(points), source(scan);
  small_gicp::estimate_covariances(target, 10); small_gicp::estimate_covariances(source, 10);
  small_gicp::GaussianVoxelMap map(.2); map.set_search_offsets(27); map.insert(target);
  small_gicp::Registration<small_gicp::GICPFactor, small_gicp::SerialReduction> registration;
  registration.optimizer.max_iterations = 50;
  const auto result = registration.align(map, source, map);
  REQUIRE(result.converged); CHECK(result.num_inliers == source.size());
  for (size_t row = 0; row < 4; ++row)
    for (size_t col = 0; col < 4; ++col) CHECK(std::abs(result.T_target_source[row][col] - expected[row][col]) < 1e-3);
}

TEST_CASE("Gaussian map retains center-cell-only search by default") {
  Cloud cloud(std::vector<std::array<double, 3>>{{{1.1, .2, .3}}});
  cloud.cov(0) = earth::mantle::math::identityMatrix<double, 4>();
  small_gicp::GaussianVoxelMap map(1.0); map.insert(cloud);
  const auto query = Point::from({.9, .2, .3, 1}); size_t index = 0; double distance = 0;
  CHECK(map.nearest_neighbor_search(query, &index, &distance) == 0);
  map.set_search_offsets(27); REQUIRE(map.nearest_neighbor_search(query, &index, &distance) == 1);
  CHECK(distance == doctest::Approx(.04));
}
