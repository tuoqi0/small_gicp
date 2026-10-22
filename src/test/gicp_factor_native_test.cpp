#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <small_gicp/util/lie.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <Eigen/LU>
#include <random>
#include <type_traits>

static_assert(std::is_same<decltype(small_gicp::GICPFactor::mahalanobis),
    earth::mantle::collection::Matrix<double,4,4>>::value,"factor precision must be native");

namespace {
struct Accept {
  template <typename... Args> bool operator()(const Args&...) const { return false; }
};
}

TEST_CASE("GICP precision and linearization agree with the original inverse") {
  std::mt19937 generator(714);
  std::normal_distribution<double> random(0, 1);
  for (int sample = 0; sample < 1000; ++sample) {
    Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
    transform.linear() = Eigen::AngleAxisd(0.37 * sample,
        Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
    transform.translation() = Eigen::Vector3d(0.1, -0.2, 0.3);
    const Eigen::Vector3d point(random(generator), random(generator), random(generator));
    small_gicp::PointCloud source(std::vector<Eigen::Vector3d>{point});
    small_gicp::PointCloud target(std::vector<Eigen::Vector3d>{
        transform * point + Eigen::Vector3d(0.01, -0.02, 0.03)});
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.23 * sample,
        Eigen::Vector3d(3, 1, -2).normalized()).toRotationMatrix();
    const double scale = std::pow(10.0, sample % 17 - 8);
    source.cov(0).setZero();
    target.cov(0).setZero();
    source.cov(0).topLeftCorner<3, 3>() =
        scale * rotation * Eigen::Vector3d(0.001, 0.1, 1).asDiagonal() * rotation.transpose();
    target.cov(0).topLeftCorner<3, 3>() = scale * Eigen::Vector3d(0.002, 0.2, 2).asDiagonal();
    small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
    small_gicp::GICPFactor factor;
    small_gicp::factor_detail::Hessian native_hessian{};
    auto hessian = small_gicp::eigen_view(native_hessian);
    small_gicp::factor_detail::Gradient native_gradient{};
    auto gradient = small_gicp::eigen_view(native_gradient);
    double error;
    REQUIRE(factor.linearize(target, source, tree, transform, 0, Accept{},
                             &native_hessian, &native_gradient, &error));
    REQUIRE(factor.inlier());
    const Eigen::Matrix4d covariance =
        target.cov(0) + transform.matrix() * source.cov(0) * transform.matrix().transpose();
    Eigen::Matrix4d inverse = Eigen::Matrix4d::Zero();
    inverse.topLeftCorner<3, 3>() = covariance.topLeftCorner<3, 3>().inverse();
    const Eigen::Map<const Eigen::Matrix<double,4,4,Eigen::RowMajor>> precision(
        factor.mahalanobis.raw_pointer());
    CHECK((precision - inverse).norm() <= 1e-11 * inverse.norm());
    Eigen::Matrix<double, 4, 6> jacobian = Eigen::Matrix<double, 4, 6>::Zero();
    jacobian.topLeftCorner<3, 3>() = transform.linear() * small_gicp::skew(point);
    jacobian.topRightCorner<3, 3>() = -transform.linear();
    const Eigen::Vector4d residual = target.point(0) - transform * source.point(0);
    const Eigen::Matrix<double, 6, 6> expected_hessian = jacobian.transpose() * inverse * jacobian;
    const Eigen::Matrix<double, 6, 1> expected_gradient = jacobian.transpose() * inverse * residual;
    const double expected_error = 0.5 * residual.dot(inverse * residual);
    CHECK((hessian - expected_hessian).norm() <= 1e-11 * expected_hessian.norm());
    CHECK((gradient - expected_gradient).norm() <= 1e-11 * expected_gradient.norm());
    CHECK(error == doctest::Approx(expected_error).epsilon(1e-11));
    // Trial poses reuse the linearization precision; they must not rebuild it.
    auto candidate=transform;
    candidate.translation()+=Eigen::Vector3d(.03,-.01,.02);
    const Eigen::Vector4d candidate_residual=target.point(0)-candidate*source.point(0);
    const double candidate_error=.5*candidate_residual.dot(inverse*candidate_residual);
    CHECK(factor.error(target,source,candidate)==doctest::Approx(candidate_error).epsilon(1e-11));
  }
}

TEST_CASE("Invalid point covariance rejects the factor without publishing partial equations") {
  small_gicp::PointCloud target(std::vector<Eigen::Vector3d>{Eigen::Vector3d::Zero()});
  small_gicp::PointCloud source(std::vector<Eigen::Vector3d>{Eigen::Vector3d::Zero()});
  small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
  small_gicp::GICPFactor factor;
  for (int invalid = 0; invalid < 3; ++invalid) {
    source.cov(0).setZero();
    target.cov(0).setZero();
    if (invalid == 1) target.cov(0)(0, 0) = -1;
    if (invalid == 2) target.cov(0)(0, 0) = std::numeric_limits<double>::quiet_NaN();
    small_gicp::factor_detail::Hessian native_hessian = small_gicp::native_matrix(Eigen::Matrix<double, 6, 6>::Constant(7));
    auto hessian = small_gicp::eigen_view(native_hessian);
    small_gicp::factor_detail::Gradient native_gradient = small_gicp::native_matrix(Eigen::Matrix<double, 6, 1>::Constant(8));
    auto gradient = small_gicp::eigen_view(native_gradient);
    double error = 9;
    CHECK_FALSE(factor.linearize(target, source, tree, Eigen::Isometry3d::Identity(), 0,
                                 Accept{}, &native_hessian, &native_gradient, &error));
    CHECK_FALSE(factor.inlier());
    CHECK((hessian.array() == 7).all());
    CHECK((gradient.array() == 8).all());
    CHECK(error == 9);
    CHECK(factor.error(target, source, Eigen::Isometry3d::Identity()) == 0);
  }
}
