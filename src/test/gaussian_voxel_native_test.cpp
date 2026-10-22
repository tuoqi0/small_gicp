#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <small_gicp/ann/gaussian_voxelmap.hpp>
#include <small_gicp/ann/sequential_voxelmap_accessor.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#include <cmath>
#include <random>
#include <type_traits>

namespace {
using Voxel = small_gicp::GaussianVoxel;
static_assert(std::is_same<decltype(Voxel::mean),
    earth::mantle::collection::Matrix<double, 4, 1>>::value,
    "voxel mean must own native storage");
static_assert(std::is_same<decltype(Voxel::cov),
    earth::mantle::collection::Matrix<double, 4, 4>>::value,
    "voxel covariance must own native storage");

// Independent reference retains the original Eigen sum/average operations.
struct ReferenceVoxel {
  bool finalized{false};
  size_t count{0};
  Eigen::Vector4d mean{Eigen::Vector4d::Zero()};
  Eigen::Matrix4d covariance{Eigen::Matrix4d::Zero()};

  void add(const Eigen::Vector4d& point, const Eigen::Matrix4d& cov,
           const Eigen::Isometry3d& transform) {
    if (finalized) {
      finalized = false;
      mean *= count;
      covariance *= count;
    }
    ++count;
    mean += point;
    covariance += transform.matrix() * cov * transform.matrix().transpose();
  }

  void finalize() {
    if (finalized) return;
    mean /= count;
    covariance /= count;
    finalized = true;
  }
};

void checkVoxel(const Voxel& voxel, const ReferenceVoxel& reference) {
  CHECK(voxel.finalized == reference.finalized);
  CHECK(voxel.num_points == reference.count);
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(voxel, 0)).isApprox(reference.mean, 1e-12));
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::cov(voxel, 0)).isApprox(reference.covariance, 1e-12));
}

struct AcceptCorrespondence {
  template <typename... Args>
  bool operator()(const Args&...) const { return false; }
};
}

TEST_CASE("native Gaussian voxel matches original transformed accumulation") {
  std::mt19937 generator(918);
  std::normal_distribution<double> random(0, 1);
  small_gicp::PointCloud cloud;
  cloud.resize(1);
  Voxel voxel;
  ReferenceVoxel reference;
  for (size_t sample = 0; sample < 1000; ++sample) {
    const Eigen::Vector3d point(random(generator), random(generator), random(generator));
    Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
    transform.linear() = Eigen::AngleAxisd(0.17 * sample,
        Eigen::Vector3d(1, -2, 3).normalized()).toRotationMatrix();
    transform.translation() = Eigen::Vector3d(0.3, -0.7, 1.2);
    small_gicp::eigen_view(cloud.point(0)) << point, 1.0;
    small_gicp::eigen_view(cloud.cov(0)).setZero();
    small_gicp::eigen_view(cloud.cov(0)).topLeftCorner<3, 3>() =
        std::pow(10.0, static_cast<int>(sample % 9) - 4) *
        (Eigen::Vector3d(0.01, 0.1, 1).asDiagonal().toDenseMatrix());
    const Eigen::Vector4d transformed = transform * small_gicp::eigen_view(cloud.point(0));
    voxel.add(Voxel::Setting{}, small_gicp::native_matrix(transformed), cloud, 0, small_gicp::native_transform(transform));
    reference.add(transformed, small_gicp::eigen_view(cloud.cov(0)), transform);
    checkVoxel(voxel, reference);
    if (sample % 7 == 0) {
      voxel.finalize();
      reference.finalize();
      checkVoxel(voxel, reference);
      voxel.finalize();
      checkVoxel(voxel, reference);
    }
  }
  voxel.finalize();
  reference.finalize();
  checkVoxel(voxel, reference);
}

TEST_CASE("native Gaussian voxel preserves covariance layout and value snapshots") {
  small_gicp::PointCloud cloud;
  cloud.resize(1);
  small_gicp::eigen_view(cloud.point(0)) = Eigen::Vector4d(1, 2, 3, 1);
  for (int row = 0; row < 4; ++row)
    for (int col = 0; col < 4; ++col) small_gicp::eigen_view(cloud.cov(0))(row, col) = 10 * row + col;
  Voxel voxel;
  voxel.add(Voxel::Setting{}, cloud.point(0), cloud, 0, small_gicp::native_transform(Eigen::Isometry3d::Identity()));
  voxel.finalize();
  CHECK(voxel.cov[0][1] == 1);
  CHECK(voxel.cov[1][0] == 10);
  const auto point_snapshot = small_gicp::eigen_matrix(small_gicp::traits::point(voxel, 0));
  const auto covariance_snapshot = small_gicp::eigen_matrix(small_gicp::traits::cov(voxel, 0));
  Voxel copy = voxel;
  copy.mean[0][0] = 9;
  CHECK(voxel.mean[0][0] == 1);
  small_gicp::eigen_view(cloud.point(0)) = Eigen::Vector4d(3, 4, 5, 1);
  voxel.add(Voxel::Setting{}, cloud.point(0), cloud, 0, small_gicp::native_transform(Eigen::Isometry3d::Identity()));
  voxel.finalize();
  CHECK(point_snapshot == Eigen::Vector4d(1, 2, 3, 1));
  CHECK(covariance_snapshot == small_gicp::eigen_view(cloud.cov(0)));
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(voxel, 0)) == Eigen::Vector4d(2, 3, 4, 1));
}

TEST_CASE("native Gaussian map supports repeated insertion search and GICP linearization") {
  small_gicp::PointCloud cloud(std::vector<Eigen::Vector3d>{{0.1, 0.2, 0.3}});
  small_gicp::eigen_view(cloud.cov(0)).setIdentity();
  small_gicp::eigen_view(cloud.cov(0))(3, 3) = 0;
  small_gicp::GaussianVoxelMap map(1.0);
  map.insert(cloud);
  small_gicp::eigen_view(cloud.point(0))[0] = 0.3;
  map.insert(cloud);
  REQUIRE(map.size() == 1);
  const small_gicp::SequentialVoxelMapAccessor<small_gicp::GaussianVoxelMap> accessor(map);
  REQUIRE(small_gicp::traits::size(accessor) == 1);
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(accessor, 0)).isApprox(Eigen::Vector4d(0.2, 0.2, 0.3, 1)));
  size_t index = 0;
  double distance = 0;
  REQUIRE(map.nearest_neighbor_search(small_gicp::eigen_view(cloud.point(0)), &index, &distance) == 1);
  CHECK(distance == doctest::Approx(0.01));
  small_gicp::GICPFactor factor;
  small_gicp::factor_detail::Hessian native_hessian{};
  auto hessian = small_gicp::eigen_view(native_hessian);
  small_gicp::factor_detail::Gradient native_gradient{};
  auto gradient = small_gicp::eigen_view(native_gradient);
  double error = 0;
  REQUIRE(factor.linearize(map, cloud, map, small_gicp::native_transform(Eigen::Isometry3d::Identity()), 0,
      AcceptCorrespondence{}, &native_hessian, &native_gradient, &error));
  CHECK(error == doctest::Approx(0.0025));
  CHECK(factor.error(map, cloud, small_gicp::native_transform(Eigen::Isometry3d::Identity())) == doctest::Approx(error));
}

TEST_CASE("native Gaussian map registration recovers a rigid pose") {
  std::vector<Eigen::Vector3d> points, scan;
  for (int row = -6; row <= 6; ++row) {
    for (int col = -6; col <= 6; ++col) {
      points.emplace_back(4, row * 0.25, col * 0.25);
      points.emplace_back(row * 0.25, 3, col * 0.25);
      points.emplace_back(row * 0.25, col * 0.25, -1);
    }
  }
  Eigen::Isometry3d expected = Eigen::Isometry3d::Identity();
  expected.translation() = Eigen::Vector3d(0.08, -0.06, 0.04);
  expected.linear() = Eigen::AngleAxisd(0.015,
      Eigen::Vector3d(1, 2, -1).normalized()).toRotationMatrix();
  for (const auto& point : points) scan.push_back(expected.inverse() * point);
  small_gicp::PointCloud target(points), source(scan);
  small_gicp::estimate_covariances(target, 10);
  small_gicp::estimate_covariances(source, 10);
  small_gicp::GaussianVoxelMap map(0.2);
  // This fixture crosses voxel boundaries. Explicitly enable the existing
  // neighborhood search; production's center-cell-only default is unchanged.
  map.set_search_offsets(27);
  map.insert(target);
  small_gicp::Registration<small_gicp::GICPFactor, small_gicp::SerialReduction> registration;
  registration.optimizer.max_iterations = 50;
  const auto result = registration.align(map, source, map);
  INFO("inliers=", result.num_inliers, ", source=", source.size(),
       ", translation=", small_gicp::eigen_transform(result.T_target_source).translation().transpose(),
       ", expected=", expected.translation().transpose());
  REQUIRE(result.converged);
  CHECK(result.num_inliers == source.size());
  CHECK((small_gicp::eigen_transform(result.T_target_source).translation() - expected.translation()).norm() < 1e-3);
  CHECK((small_gicp::eigen_transform(result.T_target_source).linear() - expected.linear()).norm() < 1e-3);
}

TEST_CASE("Gaussian map retains center-cell-only search by default") {
  small_gicp::PointCloud cloud(std::vector<Eigen::Vector3d>{{1.1, 0.2, 0.3}});
  small_gicp::eigen_view(cloud.cov(0)).setIdentity();
  small_gicp::GaussianVoxelMap map(1.0);
  map.insert(cloud);
  const Eigen::Vector4d query(0.9, 0.2, 0.3, 1);
  size_t index = 0;
  double distance = 0;
  CHECK(map.nearest_neighbor_search(query, &index, &distance) == 0);
  map.set_search_offsets(27);
  REQUIRE(map.nearest_neighbor_search(query, &index, &distance) == 1);
  CHECK(distance == doctest::Approx(0.04));
}
