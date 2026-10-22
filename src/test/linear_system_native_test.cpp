#include <robot_core/test/doctest.h>
#include <small_gicp/registration/linear_system.hpp>
#include <small_gicp/registration/termination_criteria.hpp>
#include <limits>

TEST_CASE("native damped equations preserve lower triangle sign and failure output") {
  using namespace small_gicp;
  for (int sample = 0; sample < 1000; ++sample) {
    factor_detail::Hessian hessian{};
    factor_detail::Gradient expected{}, gradient{}, delta{};
    for (size_t row = 0; row < 6; ++row) {
      expected[row][0] = std::sin(.17 * (sample + 1) * (row + 1));
      for (size_t col = 0; col <= row; ++col)
        hessian[row][col] = row == col ? 8 + row : .01 * (row + col + 1);
      for (size_t col = row + 1; col < 6; ++col) hessian[row][col] = 99;
    }
    const double damping = sample % 3 == 0 ? 0 : 1e-3;
    for (size_t row = 0; row < 6; ++row) {
      for (size_t col = 0; col < 6; ++col)
        gradient[row][0] -= hessian[std::max(row, col)][std::min(row, col)] * expected[col][0];
      gradient[row][0] -= damping * expected[row][0];
    }
    REQUIRE(detail::solve_damped_system(hessian, gradient, damping, delta));
    for (size_t row = 0; row < 6; ++row) CHECK(delta[row][0] == doctest::Approx(expected[row][0]).epsilon(1e-12));
  }
  factor_detail::Hessian invalid{};
  factor_detail::Gradient gradient{}, delta{};
  delta[0][0] = 7;
  CHECK_FALSE(detail::solve_damped_system(invalid, gradient, 0, delta));
  CHECK(delta[0][0] == 7);
  invalid[0][0] = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(detail::solve_damped_system(invalid, gradient, 1, delta));
  CHECK(delta[0][0] == 7);
}

TEST_CASE("native convergence keeps inclusive independent rotation and translation limits") {
  using namespace small_gicp;
  const TerminationCriteria criteria;
  factor_detail::Gradient delta{};
  CHECK(criteria.converged(delta));
  delta[0][0] = criteria.rotation_eps;
  delta[3][0] = criteria.translation_eps;
  CHECK(criteria.converged(delta));
  delta[0][0] *= 1.01;
  CHECK_FALSE(criteria.converged(delta));
  delta[0][0] = 0;
  delta[3][0] *= 1.01;
  CHECK_FALSE(criteria.converged(delta));
  delta[3][0] = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(criteria.converged(delta));
}
