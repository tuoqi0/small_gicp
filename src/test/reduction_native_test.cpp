#include <robot_core/test/doctest.h>
#include <small_gicp/util/eigen.hpp>
#include <small_gicp/registration/reduction.hpp>
#include <small_gicp/registration/reduction_omp.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#ifdef BUILD_WITH_TBB
#include <small_gicp/registration/reduction_tbb.hpp>
#include <tbb/global_control.h>
#endif

namespace {
struct SyntheticFactor {
  bool accepted{true};
  double weight{0};

  template <typename Rejector>
  bool linearize(const int&, const int&, const int&, const Eigen::Isometry3d&,
                 size_t, const Rejector&, small_gicp::factor_detail::Hessian* hessian,
                 small_gicp::factor_detail::Gradient* gradient, double* error) {
    // Rejected factors deliberately write outputs: reducers must skip them.
    for (int row = 0; row < 6; ++row) {
      (*gradient)[row][0] = weight * (row - 3);
      for (int col = 0; col < 6; ++col)
        (*hessian)[row][col] = weight * (10 * row + col + 1);
    }
    *error = 0.5 * weight;
    return accepted;
  }

  double error(const int&, const int&, const Eigen::Isometry3d&) const {
    return accepted ? 0.5 * weight : 0;
  }
};

template <typename Reduction>
void checkReduction(const Reduction& reduction, size_t count, bool reject_all) {
  std::vector<SyntheticFactor> factors;
  Eigen::Matrix<double, 6, 6> expected_hessian = Eigen::Matrix<double, 6, 6>::Zero();
  Eigen::Matrix<double, 6, 1> expected_gradient = Eigen::Matrix<double, 6, 1>::Zero();
  double expected_error = 0;
  for (size_t i = 0; i < count; ++i) {
    // Binary-exact contributions allow strict equality regardless of scheduling.
    const double weight = static_cast<double>(i + 1) * 0.125;
    const bool accepted = !reject_all && i % 3 != 0;
    factors.push_back(SyntheticFactor{accepted, weight});
    if (!accepted) continue;
    small_gicp::factor_detail::Hessian hessian{};
    small_gicp::factor_detail::Gradient gradient{};
    double error = 0;
    factors.back().linearize(0, 0, 0, Eigen::Isometry3d::Identity(), i, 0,
                            &hessian, &gradient, &error);
    expected_hessian += small_gicp::eigen_view(hessian);
    expected_gradient += small_gicp::eigen_view(gradient);
    expected_error += error;
  }
  const auto result = reduction.linearize(0, 0, 0, 0, Eigen::Isometry3d::Identity(), factors);
  CHECK(small_gicp::eigen_view(std::get<0>(result)) == expected_hessian);
  CHECK(small_gicp::eigen_view(std::get<1>(result)) == expected_gradient);
  CHECK(std::get<2>(result) == expected_error);
  CHECK(reduction.error(0, 0, Eigen::Isometry3d::Identity(), factors) == expected_error);
}
}

TEST_CASE("native serial reduction preserves accepted rejected and empty sums") {
  const small_gicp::SerialReduction reduction;
  checkReduction(reduction, 0, false);
  checkReduction(reduction, 97, false);
  checkReduction(reduction, 97, true);
}

TEST_CASE("native OpenMP reduction preserves layout and sums across thread counts") {
  small_gicp::ParallelReductionOMP reduction;
  for (int threads : {1, 2, 4}) {
    reduction.num_threads = threads;
    checkReduction(reduction, 0, false);
    checkReduction(reduction, 97, false);
    checkReduction(reduction, 97, true);
  }
}

#ifdef BUILD_WITH_TBB
TEST_CASE("native TBB reduction preserves accepted rejected and empty sums") {
  const tbb::global_control limit(tbb::global_control::max_allowed_parallelism, 2);
  const small_gicp::ParallelReductionTBB reduction;
  checkReduction(reduction, 0, false);
  checkReduction(reduction, 97, false);
  checkReduction(reduction, 97, true);
}
#endif

TEST_CASE("native OpenMP reduction retains GICP pose recovery") {
  std::vector<Eigen::Vector3d> map_points, scan_points;
  const Eigen::Vector3d offset(0.15, -0.1, 0.05);
  for (int row = -6; row <= 6; ++row) {
    for (int col = -6; col <= 6; ++col) {
      map_points.emplace_back(4, row * 0.2, col * 0.2);
      map_points.emplace_back(row * 0.2, 3, col * 0.2);
      map_points.emplace_back(row * 0.2, col * 0.2, -1);
    }
  }
  for (const auto& point : map_points) scan_points.push_back(point - offset);
  small_gicp::PointCloud target(map_points), source(scan_points);
  small_gicp::estimate_covariances(target, 10);
  small_gicp::estimate_covariances(source, 10);
  small_gicp::UnsafeKdTree<small_gicp::PointCloud> tree(target);
  small_gicp::Registration<small_gicp::GICPFactor, small_gicp::ParallelReductionOMP> registration;
  registration.optimizer.max_iterations = 50;
  for (int threads : {1, 2, 4}) {
    registration.reduction.num_threads = threads;
    const auto result = registration.align(target, source, tree);
    REQUIRE(result.converged);
    CHECK(result.num_inliers == source.size());
    CHECK((result.T_target_source.translation() - offset).norm() < 1e-3);
    CHECK((result.T_target_source.linear() - Eigen::Matrix3d::Identity()).norm() < 1e-3);
  }
}
