#include <robot_core/test/doctest.h>
#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/registration/reduction_omp.hpp>
#include <small_gicp/factors/icp_factor.hpp>
#include <small_gicp/factors/plane_icp_factor.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#include <small_gicp/util/downsampling.hpp>
#include <small_gicp/ann/gaussian_voxelmap.hpp>
#include <small_gicp/ann/projective_search.hpp>
#include <array>
#include <vector>

using NativeCloud = small_gicp::PointCloud;

namespace {
template <typename Factor, typename Reduction, typename Optimizer = small_gicp::LevenbergMarquardtOptimizer>
void checkNativeRegistration() {
  using namespace small_gicp;
  using earth::mantle::math::operator*;
  using earth::mantle::math::operator+;
  std::vector<std::array<double, 4>> target_points, source_points;
  for (int row = -5; row <= 5; ++row) {
    for (int col = -5; col <= 5; ++col) {
      target_points.push_back({{4, row * .25, col * .25, 1}});
      target_points.push_back({{row * .25, 3, col * .25, 1}});
      target_points.push_back({{row * .25, col * .25, -1, 1}});
    }
  }
  const auto expected = se3_exp(factor_detail::Gradient::from({.03, -.02, .01, .08, -.06, .04}));
  const auto inverse_rotation = earth::mantle::math::transpose(rotation_matrix(expected));
  for (const auto& point : target_points) {
    const auto shifted = earth::mantle::collection::Matrix<double, 3, 1>::from({
        point[0] - expected[0][3], point[1] - expected[1][3], point[2] - expected[2][3]});
    const auto transformed = inverse_rotation * shifted;
    source_points.push_back({{transformed[0][0], transformed[1][0], transformed[2][0], 1}});
  }
  NativeCloud target(target_points), source(source_points);
  estimate_normals_covariances(target, 10);
  estimate_normals_covariances(source, 10);
  const UnsafeKdTree<NativeCloud> tree(target);
  Registration<Factor, Reduction, NullFactor, DistanceRejector, Optimizer> registration;
  registration.optimizer.max_iterations = 50;
  const auto result = registration.align(target, source, tree);
  REQUIRE(result.converged);
  CHECK(result.num_inliers == source.size());
  for (size_t row = 0; row < 4; ++row)
    for (size_t col = 0; col < 4; ++col)
      CHECK(std::abs(result.T_target_source[row][col] - expected[row][col]) < 1e-3);
}
}

TEST_CASE("registration core recovers rigid poses without Eigen PCL or Boost headers") {
  checkNativeRegistration<small_gicp::ICPFactor, small_gicp::SerialReduction>();
  checkNativeRegistration<small_gicp::PointToPlaneICPFactor, small_gicp::SerialReduction>();
  checkNativeRegistration<small_gicp::GICPFactor, small_gicp::SerialReduction>();
  checkNativeRegistration<small_gicp::GICPFactor, small_gicp::SerialReduction, small_gicp::GaussNewtonOptimizer>();
}

TEST_CASE("actual point cloud preprocessing and voxel lookup have no Eigen storage boundary") {
  NativeCloud cloud(std::vector<std::array<double, 3>>{{{-.15, .1, .2}}, {{-.12, .1, .2}}, {{.15, .1, .2}}});
  auto sampled = small_gicp::voxelgrid_sampling(cloud, .1);
  REQUIRE(sampled->size() == 2);
  small_gicp::estimate_covariances(*sampled, 5);
  small_gicp::GaussianVoxelMap map(.5);
  map.insert(*sampled);
  size_t index = 0;
  double distance = 0;
  REQUIRE(map.nearest_neighbor_search(std::array<double, 4>{{-.135, .1, .2, 1}}, &index, &distance) == 1);
  CHECK(distance < 1e-12);
  const auto covariance = small_gicp::traits::cov(map, index);
  CHECK(covariance[0][0] == 1);
}

TEST_CASE("native projective lookup preserves pixel overwrite and horizontal wrap") {
  const NativeCloud cloud(std::vector<std::array<double, 3>>{{{0, 0, 1}}, {{0, 0, 2}},
      {{1, 0, 0}}, {{-1e-10, 0, -1}}});
  small_gicp::UnsafeProjectiveSearch<NativeCloud> search(64, 32, cloud);
  search.search_window_h = 0;
  search.search_window_v = 0;
  size_t index = 0;
  double distance = 0;
  REQUIRE(search.nearest_neighbor_search(std::array<double, 4>{{0, 0, 1.1, 1}}, &index, &distance) == 1);
  CHECK(index == 1);  // The last point in the pixel replaces the earlier one.
  CHECK(distance == doctest::Approx(.81));
  REQUIRE(search.nearest_neighbor_search(std::array<double, 4>{{0, 0, -1, 1}}, &index, &distance) == 1);
  CHECK(index == 3);
  CHECK(distance < 1e-18);
}
