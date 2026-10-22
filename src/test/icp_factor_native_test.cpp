#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <small_gicp/util/lie.hpp>
#include <small_gicp/factors/icp_factor.hpp>
#include <small_gicp/factors/plane_icp_factor.hpp>
#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/registration/registration.hpp>
#include <random>

namespace {
struct Accept {
  template <typename... Args> bool operator()(const Args&...) const { return false; }
};
struct Reject {
  template <typename... Args> bool operator()(const Args&...) const { return true; }
};
struct MissingNeighbor {
  template <typename Query>
  size_t nearest_neighbor_search(const Query&, size_t*, double*) const { return 0; }
};
}

namespace small_gicp { namespace traits {
template <> struct Traits<MissingNeighbor> {
  template <typename Query>
  static size_t nearest_neighbor_search(const MissingNeighbor& tree,
      const Query& point, size_t* index, double* distance) {
    return tree.nearest_neighbor_search(point, index, distance);
  }
};
}}

namespace {
template <typename Factor>
void checkRejectedCorrespondences() {
  small_gicp::PointCloud source(std::vector<Eigen::Vector3d>{{0.1, 0.2, 0.3}});
  small_gicp::PointCloud target(std::vector<Eigen::Vector3d>{{0.2, 0.3, 0.4}});
  target.normal(0) = Eigen::Vector4d(1, 2, 3, 0);
  small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
  Factor factor;
  small_gicp::factor_detail::Hessian native_hessian{};
  auto hessian = small_gicp::eigen_view(native_hessian);
  small_gicp::factor_detail::Gradient native_gradient{};
  auto gradient = small_gicp::eigen_view(native_gradient);
  double error = 0;
  REQUIRE(factor.linearize(target, source, tree, small_gicp::native_transform(Eigen::Isometry3d::Identity()), 0,
      Accept{}, &native_hessian, &native_gradient, &error));
  REQUIRE(factor.inlier());
  hessian.setConstant(7);
  gradient.setConstant(8);
  error = 9;
  CHECK_FALSE(factor.linearize(target, source, tree, small_gicp::native_transform(Eigen::Isometry3d::Identity()), 0,
      Reject{}, &native_hessian, &native_gradient, &error));
  CHECK_FALSE(factor.inlier());
  CHECK(factor.error(target, source, small_gicp::native_transform(Eigen::Isometry3d::Identity())) == 0);
  CHECK((hessian.array() == 7).all());
  CHECK((gradient.array() == 8).all());
  CHECK(error == 9);
  CHECK_FALSE(factor.linearize(target, source, MissingNeighbor{}, small_gicp::native_transform(Eigen::Isometry3d::Identity()), 0,
      Accept{}, &native_hessian, &native_gradient, &error));
  CHECK_FALSE(factor.inlier());
  CHECK((hessian.array() == 7).all());
  CHECK((gradient.array() == 8).all());
  CHECK(error == 9);
}

template <typename Factor>
void checkPoseRecovery() {
  std::vector<Eigen::Vector3d> points, scan;
  const Eigen::Vector3d offset(0.08, -0.06, 0.04);
  for (int row = -6; row <= 6; ++row) {
    for (int col = -6; col <= 6; ++col) {
      points.emplace_back(4, row * 0.25, col * 0.25);
      points.emplace_back(row * 0.25, 3, col * 0.25);
      points.emplace_back(row * 0.25, col * 0.25, -1);
    }
  }
  for (const auto& point : points) scan.push_back(point - offset);
  small_gicp::PointCloud target(points), source(scan);
  // Non-axis-aligned weights exercise the original componentwise normal rule.
  for (size_t i = 0; i < target.size(); ++i)
    target.normal(i) = Eigen::Vector4d(0.3, -0.4, 0.5, 0);
  small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
  small_gicp::Registration<Factor, small_gicp::SerialReduction> registration;
  registration.optimizer.max_iterations = 50;
  const auto result = registration.align(target, source, tree);
  REQUIRE(result.converged);
  CHECK(result.num_inliers == source.size());
  CHECK((small_gicp::eigen_transform(result.T_target_source).translation() - offset).norm() < 1e-3);
  CHECK((small_gicp::eigen_transform(result.T_target_source).linear() - Eigen::Matrix3d::Identity()).norm() < 1e-3);
}
}

TEST_CASE("native ICP equations and trial-pose error match original Eigen formulas") {
  std::mt19937 generator(621);
  std::normal_distribution<double> random(0, 1);
  for (int sample = 0; sample < 1000; ++sample) {
    const Eigen::Vector3d point(random(generator), random(generator), random(generator));
    Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
    transform.linear() = Eigen::AngleAxisd(0.13 * sample,
        Eigen::Vector3d(1, 2, -3).normalized()).toRotationMatrix();
    transform.translation() = Eigen::Vector3d(0.1, -0.3, 0.4);
    small_gicp::PointCloud source(std::vector<Eigen::Vector3d>{point});
    small_gicp::PointCloud target(std::vector<Eigen::Vector3d>{
        transform * point + Eigen::Vector3d(0.03, -0.04, 0.05)});
    small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
    small_gicp::ICPFactor factor;
    small_gicp::factor_detail::Hessian native_hessian{};
    auto hessian = small_gicp::eigen_view(native_hessian);
    small_gicp::factor_detail::Gradient native_gradient{};
    auto gradient = small_gicp::eigen_view(native_gradient);
    double error = 0;
    REQUIRE(factor.linearize(target, source, tree, small_gicp::native_transform(transform), 0, Accept{}, &native_hessian, &native_gradient, &error));
    const Eigen::Vector4d residual = target.point(0) - transform * source.point(0);
    Eigen::Matrix<double, 4, 6> jacobian = Eigen::Matrix<double, 4, 6>::Zero();
    jacobian.topLeftCorner<3, 3>() = transform.linear() * small_gicp::skew(point);
    jacobian.topRightCorner<3, 3>() = -transform.linear();
    CHECK(hessian.isApprox(jacobian.transpose() * jacobian, 1e-12));
    CHECK(gradient.isApprox(jacobian.transpose() * residual, 1e-12));
    CHECK(error == doctest::Approx(0.5 * residual.squaredNorm()).epsilon(1e-12));
    auto trial = transform;
    trial.translation() += Eigen::Vector3d(0.3, -0.1, 0.2);
    const Eigen::Vector4d trial_residual = target.point(0) - trial * source.point(0);
    CHECK(factor.error(target, source, small_gicp::native_transform(trial)) == doctest::Approx(0.5 * trial_residual.squaredNorm()).epsilon(1e-12));
  }
}

TEST_CASE("native normal-weighted ICP preserves componentwise not scalar plane error") {
  std::mt19937 generator(914);
  std::normal_distribution<double> random(0, 1);
  for (int sample = 0; sample < 1000; ++sample) {
    const Eigen::Vector3d point(random(generator), random(generator), random(generator));
    Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
    transform.linear() = Eigen::AngleAxisd(0.19 * sample,
        Eigen::Vector3d(2, -1, 3).normalized()).toRotationMatrix();
    transform.translation() = Eigen::Vector3d(0.1, -0.3, 0.4);
    small_gicp::PointCloud source(std::vector<Eigen::Vector3d>{point});
    small_gicp::PointCloud target(std::vector<Eigen::Vector3d>{
        transform * point + Eigen::Vector3d(0.03, -0.04, 0.05)});
    target.normal(0) = Eigen::Vector4d(0.3, -0.4, 0.5, 0);
    small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
    small_gicp::PointToPlaneICPFactor factor;
    small_gicp::factor_detail::Hessian native_hessian{};
    auto hessian = small_gicp::eigen_view(native_hessian);
    small_gicp::factor_detail::Gradient native_gradient{};
    auto gradient = small_gicp::eigen_view(native_gradient);
    double error = 0;
    REQUIRE(factor.linearize(target, source, tree, small_gicp::native_transform(transform), 0, Accept{}, &native_hessian, &native_gradient, &error));
    const Eigen::Vector4d normal = target.normal(0);
    const Eigen::Vector4d residual = target.point(0) - transform * source.point(0);
    const Eigen::Vector4d weighted = normal.array() * residual.array();
    Eigen::Matrix<double, 4, 6> jacobian = Eigen::Matrix<double, 4, 6>::Zero();
    jacobian.topLeftCorner<3, 3>() = normal.head<3>().asDiagonal() * transform.linear() * small_gicp::skew(point);
    jacobian.topRightCorner<3, 3>() = normal.head<3>().asDiagonal() * (-transform.linear());
    CHECK(hessian.isApprox(jacobian.transpose() * jacobian, 1e-12));
    CHECK(gradient.isApprox(jacobian.transpose() * weighted, 1e-12));
    CHECK(error == doctest::Approx(0.5 * weighted.squaredNorm()).epsilon(1e-12));
    CHECK(std::abs(error - 0.5 * std::pow(normal.dot(residual), 2)) > 1e-4);
    auto trial = transform;
    trial.translation() += Eigen::Vector3d(0.3, -0.1, 0.2);
    const Eigen::Vector4d trial_residual = target.point(0) - trial * source.point(0);
    const Eigen::Vector4d trial_weighted = normal.array() * trial_residual.array();
    CHECK(factor.error(target, source, small_gicp::native_transform(trial)) == doctest::Approx(0.5 * trial_weighted.squaredNorm()).epsilon(1e-12));
  }
}

TEST_CASE("native ICP factors preserve correspondence rejection and untouched outputs") {
  checkRejectedCorrespondences<small_gicp::ICPFactor>();
  checkRejectedCorrespondences<small_gicp::PointToPlaneICPFactor>();
}

TEST_CASE("native ICP factors retain iterative pose recovery") {
  checkPoseRecovery<small_gicp::ICPFactor>();
  checkPoseRecovery<small_gicp::PointToPlaneICPFactor>();
}

TEST_CASE("native normal-weighted ICP retains zero-normal behavior") {
  small_gicp::PointCloud source(std::vector<Eigen::Vector3d>{{1, 2, 3}});
  small_gicp::PointCloud target(std::vector<Eigen::Vector3d>{{1.1, 2.1, 3.1}});
  target.normal(0).setZero();
  small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
  small_gicp::PointToPlaneICPFactor factor;
  small_gicp::factor_detail::Hessian native_hessian{};
  auto hessian = small_gicp::eigen_view(native_hessian);
  small_gicp::factor_detail::Gradient native_gradient{};
  auto gradient = small_gicp::eigen_view(native_gradient);
  double error = -1;
  REQUIRE(factor.linearize(target, source, tree, small_gicp::native_transform(Eigen::Isometry3d::Identity()), 0,
      Accept{}, &native_hessian, &native_gradient, &error));
  CHECK(factor.inlier());
  CHECK(hessian.isZero());
  CHECK(gradient.isZero());
  CHECK(error == 0);
  CHECK(factor.error(target, source, small_gicp::native_transform(Eigen::Isometry3d::Identity())) == 0);
}
