#include <robot_core/test/doctest.h>
#include <small_gicp/factors/linearization.hpp>
#include <limits>

TEST_CASE("native factor equations retain rotation first and homogeneous padding") {
  using namespace small_gicp::factor_detail;
  const auto identity = earth::mantle::math::identityMatrix<double, 3>();
  const auto jacobian = pointJacobian(Point::from({1, 2, 3, 1}), identity);
  const double expected[4][6] = {
      {0, -3, 2, -1, 0, 0}, {3, 0, -1, 0, -1, 0},
      {-2, 1, 0, 0, 0, -1}, {0, 0, 0, 0, 0, 0}};
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 6; ++col) CHECK(jacobian[row][col] == expected[row][col]);
  const auto residual = Point::from({0.1, -0.2, 0.3, 0});
  Hessian hessian{}, weighted_hessian{};
  Gradient gradient{}, weighted_gradient{};
  double error = 0, weighted_error = 0;
  linearizeResidual(jacobian, residual, hessian, gradient, error);
  auto precision = earth::mantle::math::identityMatrix<double, 4>();
  precision[3][3] = 0;
  linearizeResidual(jacobian, residual, precision, weighted_hessian, weighted_gradient, weighted_error);
  CHECK(error == doctest::Approx(0.07));
  CHECK(weighted_error == error);
  for (std::size_t row = 0; row < 6; ++row) {
    double expected_gradient = 0;
    for (std::size_t axis = 0; axis < 4; ++axis)
      expected_gradient += expected[axis][row] * residual[axis][0];
    CHECK(gradient[row][0] == doctest::Approx(expected_gradient));
    CHECK(weighted_gradient[row][0] == gradient[row][0]);
    for (std::size_t col = 0; col < 6; ++col) {
      double expected_hessian = 0;
      for (std::size_t axis = 0; axis < 4; ++axis)
        expected_hessian += expected[axis][row] * expected[axis][col];
      CHECK(hessian[row][col] == expected_hessian);
      CHECK(weighted_hessian[row][col] == hessian[row][col]);
    }
  }
}

TEST_CASE("native GICP precision rejects invalid covariance without changing the output") {
  using namespace small_gicp::factor_detail;
  Covariance covariance{}, output{};
  const auto transform = earth::mantle::math::identityMatrix<double, 4>();
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 4; ++col) output[row][col] = 7;
  CHECK_FALSE(gicpPrecision(covariance, covariance, transform, output));
  CHECK(output[0][0] == 7);
  for (std::size_t axis = 0; axis < 3; ++axis) covariance[axis][axis] = 0.5;
  REQUIRE(gicpPrecision(covariance, covariance, transform, output));
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 4; ++col)
      CHECK(output[row][col] == ((row == col && row < 3) ? 1 : 0));
  covariance[0][0] = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(gicpPrecision(covariance, covariance, transform, output));
  CHECK(output[0][0] == 1);
}
