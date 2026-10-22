#include <robot_core/test/doctest.h>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/ann/kdtree.hpp>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <numeric>
#include <random>

TEST_CASE("projection preserves Eigen sampling and homogeneous-coordinate normalization") {
  std::mt19937 generator(614);
  std::normal_distribution<double> random(0, 1);
  for (const int limit : {7, 128, 1000}) {
    for (int trial = 0; trial < 100; ++trial) {
      std::vector<Eigen::Vector3d> points;
      for (int i = 0; i < 257; ++i)
        points.emplace_back(4 + random(generator), -2 + 2 * random(generator), 3 * random(generator));
      small_gicp::PointCloud cloud(points);
      // Preserve sum(w), rather than silently replacing it with sample count.
      for (size_t i = 0; i < cloud.size(); ++i) cloud.point(i).w() = 1 + 0.01 * (i % 3);
      std::vector<size_t> indices(cloud.size());
      std::iota(indices.begin(), indices.end(), 0);
      std::shuffle(indices.begin(), indices.end(), generator);
      const size_t step = cloud.size() < static_cast<size_t>(limit) ? 1 : cloud.size() / limit;
      Eigen::Vector4d sum = Eigen::Vector4d::Zero();
      Eigen::Vector4d squares = Eigen::Vector4d::Zero();
      Eigen::Matrix4d cross = Eigen::Matrix4d::Zero();
      for (size_t i = 0; i < cloud.size() / step; ++i) {
        const auto& point = cloud.point(indices[step * i]);
        sum += point;
        squares += point.cwiseProduct(point);
        cross += point * point.transpose();
      }
      const Eigen::Vector4d mean = sum / sum.w();
      const Eigen::Vector4d variance = squares - mean.cwiseProduct(sum);
      const int expected_axis = variance[0] > variance[1] ?
        (variance[0] > variance[2] ? 0 : 2) : (variance[1] > variance[2] ? 1 : 2);
      const Eigen::Matrix4d covariance = (cross - mean * sum.transpose()) / sum.w();
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> oracle(covariance.topLeftCorner<3, 3>());
      REQUIRE(oracle.info() == Eigen::Success);
      small_gicp::ProjectionSetting setting;
      setting.max_scan_count = limit;
      const auto axis = small_gicp::AxisAlignedProjection::find_axis(
        cloud, indices.begin(), indices.end(), setting);
      const auto normal = small_gicp::NormalProjection::find_axis(
        cloud, indices.begin(), indices.end(), setting);
      CHECK(axis.axis == expected_axis);
      const Eigen::Vector3d direction = Eigen::Map<const Eigen::Vector3d>(normal.normal.data());
      CHECK(std::abs(direction.dot(oracle.eigenvectors().col(2))) > 1 - 1e-10);
      const Eigen::Vector4d query(0.4, -0.7, 1.2, 1);
      CHECK(normal(query) == doctest::Approx(direction.dot(query.head<3>())).epsilon(1e-12));
    }
  }
}

TEST_CASE("Normal projection selects the largest covariance eigenspace") {
  std::mt19937 generator(173);
  std::normal_distribution<double> random(0, 1);
  for (int sample = 0; sample < 500; ++sample) {
    const double scale = std::pow(10.0, sample % 13 - 6);
    std::vector<Eigen::Vector3d> points;
    for (int index = 0; index < 100; ++index)
      points.emplace_back(scale * random(generator), 2 * scale * random(generator),
                          3 * scale * random(generator));
    small_gicp::PointCloud cloud(points);
    std::vector<std::size_t> indices(points.size());
    std::iota(indices.begin(), indices.end(), 0);
    const auto projection = small_gicp::NormalProjection::find_axis(
        cloud, indices.begin(), indices.end(), small_gicp::ProjectionSetting{});
    const Eigen::Vector3d direction = Eigen::Map<const Eigen::Vector3d>(projection.normal.data());
    Eigen::Vector3d sum = Eigen::Vector3d::Zero();
    Eigen::Matrix3d square = Eigen::Matrix3d::Zero();
    for (const auto& point : points) { sum += point; square += point * point.transpose(); }
    const Eigen::Matrix3d covariance =
        (square - sum * sum.transpose() / static_cast<double>(points.size())) /
        static_cast<double>(points.size());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> oracle(covariance);
    REQUIRE(oracle.info() == Eigen::Success);
    CHECK(direction.norm() == doctest::Approx(1).epsilon(1e-12));
    CHECK(std::abs(direction.dot(oracle.eigenvectors().col(2))) > 1 - 1e-10);
    CHECK((covariance * direction - oracle.eigenvalues()[2] * direction).norm()
          <= 1e-11 * covariance.norm());
  }
}

TEST_CASE("Normal projection KD-tree neighbors agree with exhaustive search") {
  std::vector<Eigen::Vector3d> points;
  for (int index = 0; index < 257; ++index)
    points.emplace_back(std::sin(0.13 * index), std::cos(0.23 * index), 0.01 * index);
  small_gicp::PointCloud cloud(points);
  small_gicp::UnsafeKdTree<small_gicp::PointCloud, small_gicp::NormalProjection> tree(cloud);
  for (int index = 0; index < 100; ++index) {
    const Eigen::Vector4d query(0.1 * index, std::sin(0.17 * index), -0.03 * index, 1);
    std::vector<double> expected;
    for (const auto& point : points) expected.push_back((point - query.head<3>()).squaredNorm());
    std::sort(expected.begin(), expected.end());
    std::size_t neighbors[5];
    double distances[5];
    REQUIRE(tree.knn_search(query, 5, neighbors, distances) == 5);
    for (std::size_t k = 0; k < 5; ++k)
      CHECK(distances[k] == doctest::Approx(expected[k]).epsilon(1e-12));
  }
}

TEST_CASE("Duplicate points give a finite unit split direction") {
  small_gicp::PointCloud cloud(std::vector<Eigen::Vector3d>(64, Eigen::Vector3d(2, 3, 4)));
  std::vector<std::size_t> indices(64);
  std::iota(indices.begin(), indices.end(), 0);
  const auto projection = small_gicp::NormalProjection::find_axis(
      cloud, indices.begin(), indices.end(), small_gicp::ProjectionSetting{});
  CHECK(Eigen::Map<const Eigen::Vector3d>(projection.normal.data()).norm() == doctest::Approx(1));
}
