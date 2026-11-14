#include <cmath>
#include <tuple>
#include <vector>
#include <robot_core/test/doctest.h>
#include <small_gicp/registration/optimizer.hpp>
#include <small_gicp/registration/termination_criteria.hpp>
#include <small_gicp/factors/general_factor.hpp>

namespace {
struct Inlier {
  bool inlier() const { return true; }
};
struct QuadraticReduction {
  double diagonal{1};
  double error(const int&, const int&, const small_gicp::Transform& pose, std::vector<Inlier>&) const {
    const double dx = pose[0][3] - .2;
    return .5 * dx * dx;
  }
  template <typename Rejector>
  auto linearize(const int& target, const int& source, const int&, const Rejector&,
      const small_gicp::Transform& pose, std::vector<Inlier>& factors) const {
    small_gicp::factor_detail::Hessian hessian{};
    small_gicp::factor_detail::Gradient gradient{};
    for (size_t axis = 0; axis < 6; ++axis) hessian[axis][axis] = diagonal;
    gradient[3][0] = pose[0][3] - .2;
    return std::make_tuple(hessian, gradient, error(target, source, pose, factors));
  }
};

template <typename Optimizer>
void checkOptimizer(Optimizer optimizer) {
  QuadraticReduction reduction;
  small_gicp::NullFactor general;
  std::vector<Inlier> factors(2);
  const auto result = optimizer.optimize(0, 0, 0, 0, small_gicp::TerminationCriteria{},
      reduction, earth::mantle::math::identityMatrix<double, 4>(), factors, general);
  CHECK(result.converged);
  CHECK(result.num_inliers == 2);
  CHECK(result.T_target_source[0][3] == doctest::Approx(.2).epsilon(1e-6));
  CHECK(result.H[3][3] == 1);
  CHECK(std::isfinite(result.b[3][0]));
}
}

TEST_CASE("GN and LM consume native factor systems and keep pose update semantics") {
  checkOptimizer(small_gicp::GaussNewtonOptimizer{});
  checkOptimizer(small_gicp::LevenbergMarquardtOptimizer{});
}

TEST_CASE("native optimizer failure keeps initial pose and reports the attempted system") {
  QuadraticReduction reduction{0};
  small_gicp::NullFactor general;
  std::vector<Inlier> factors(2);
  small_gicp::GaussNewtonOptimizer gn;
  gn.lambda = 0;
  const auto initial = earth::mantle::math::identityMatrix<double, 4>();
  const auto failed = gn.optimize(0, 0, 0, 0, small_gicp::TerminationCriteria{}, reduction, initial, factors, general);
  CHECK_FALSE(failed.converged);
  for(size_t row=0;row<4;++row)
    for(size_t col=0;col<4;++col) CHECK(failed.T_target_source[row][col]==initial[row][col]);
  CHECK(failed.iterations == 0);
  CHECK(failed.H[3][3] == 0);
  CHECK(failed.b[3][0] == -.2);
  CHECK(failed.num_inliers == 2);

  reduction.diagonal = -1;
  small_gicp::LevenbergMarquardtOptimizer lm;
  lm.max_inner_iterations = 3;
  const auto rejected = lm.optimize(0, 0, 0, 0, small_gicp::TerminationCriteria{}, reduction, initial, factors, general);
  CHECK_FALSE(rejected.converged);
  for(size_t row=0;row<4;++row)
    for(size_t col=0;col<4;++col) CHECK(rejected.T_target_source[row][col]==initial[row][col]);
  CHECK(rejected.iterations == 0);
  CHECK(rejected.H[3][3] == -1);
}

TEST_CASE("native DoF masks retain soft diagonal penalties without modifying gradient") {
  small_gicp::RestrictDoFFactor factor;
  factor.set_rotation_mask(std::array<double, 3>{{0, 1, 0}});
  factor.set_translation_mask(std::array<double, 3>{{1, 0, 1}});
  auto hessian = earth::mantle::math::identityMatrix<double, 6>();
  small_gicp::factor_detail::Gradient gradient{};
  gradient[0][0] = 3;
  double error = .5;
  factor.update_linearized_system(0, 0, 0, earth::mantle::math::identityMatrix<double, 4>(), &hessian, &gradient, &error);
  for (size_t axis = 0; axis < 6; ++axis)
    CHECK(hessian[axis][axis] == (axis == 0 || axis == 2 || axis == 4 ? 1 + factor.lambda : 1));
  CHECK(gradient[0][0] == 3);
  CHECK(error == .5);
}
