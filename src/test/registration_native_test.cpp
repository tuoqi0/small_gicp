#include <robot_core/test/doctest.h>
#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/registration/reduction_omp.hpp>
#include <small_gicp/factors/icp_factor.hpp>
#include <small_gicp/factors/plane_icp_factor.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <array>
#include <vector>

using NativeCloud = std::vector<std::array<double, 4>>;
namespace small_gicp { namespace traits {
template <> struct Traits<NativeCloud> {
  static size_t size(const NativeCloud& points) { return points.size(); }
  static std::array<double, 4> point(const NativeCloud& points, size_t i) { return points[i]; }
  static std::array<double, 4> normal(const NativeCloud&, size_t) { return {{.3, -.4, .5, 0}}; }
  static factor_detail::Covariance cov(const NativeCloud&, size_t) {
    factor_detail::Covariance covariance{};
    for (size_t axis = 0; axis < 3; ++axis) covariance[axis][axis] = 1;
    return covariance;
  }
};
}}

namespace {
template <typename Factor, typename Reduction, typename Optimizer = small_gicp::LevenbergMarquardtOptimizer>
void checkNativeRegistration() {
  using namespace small_gicp;
  using earth::mantle::math::operator*;
  using earth::mantle::math::operator+;
  NativeCloud target, source;
  for (int row = -5; row <= 5; ++row) {
    for (int col = -5; col <= 5; ++col) {
      target.push_back({{4, row * .25, col * .25, 1}});
      target.push_back({{row * .25, 3, col * .25, 1}});
      target.push_back({{row * .25, col * .25, -1, 1}});
    }
  }
  const auto expected = se3_exp(factor_detail::Gradient::from({.03, -.02, .01, .08, -.06, .04}));
  const auto inverse_rotation = earth::mantle::math::transpose(rotation_matrix(expected));
  for (const auto& point : target) {
    const auto shifted = earth::mantle::collection::Matrix<double, 3, 1>::from({
        point[0] - expected[0][3], point[1] - expected[1][3], point[2] - expected[2][3]});
    const auto transformed = inverse_rotation * shifted;
    source.push_back({{transformed[0][0], transformed[1][0], transformed[2][0], 1}});
  }
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
