#include <robot_core/test/doctest.h>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#include <random>
#include <limits>

namespace {
struct OrderedNeighbors {
  size_t count;
  size_t knn_search(const Eigen::Vector4d&, int, size_t* indices, double* distances) const {
    for (size_t i = 0; i < count; ++i) {
      indices[i] = i;
      distances[i] = 0;
    }
    return count;
  }
};
}

TEST_CASE("native local features preserve the Eigen population covariance formula") {
  std::mt19937 generator(412);
  std::normal_distribution<double> sample(0, 1);
  for (const size_t count : {5U, 20U, 64U}) {
    for (int trial = 0; trial < 100; ++trial) {
      std::vector<Eigen::Vector3d> points;
      for (size_t i = 0; i < count; ++i)
        points.emplace_back(4 + sample(generator), -2 + sample(generator), 3 + 0.1 * sample(generator));
      small_gicp::PointCloud cloud(points);
      Eigen::Vector4d sum = Eigen::Vector4d::Zero();
      Eigen::Matrix4d cross = Eigen::Matrix4d::Zero();
      for (size_t i = 0; i < count; ++i) {
        sum += cloud.point(i);
        cross += cloud.point(i) * cloud.point(i).transpose();
      }
      const Eigen::Vector4d mean = sum / count;
      const Eigen::Matrix4d covariance = (cross - mean * sum.transpose()) / count;
      std::array<double, 9> lower{};
      for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
          lower[row * 3 + col] = covariance(std::max(row, col), std::min(row, col));
      earth::mantle::math::SymmetricEigen3 expected;
      REQUIRE(earth::mantle::math::symmetricEigen3(lower, expected));
      Eigen::Matrix3d vectors;
      for (int row = 0; row < 3; ++row)
        for (int col = 0; col < 3; ++col)
          vectors(row, col) = expected.vectors[row * 3 + col];
      Eigen::Vector4d normal = Eigen::Vector4d::Zero();
      normal.head<3>() = vectors.col(0).normalized();
      if (cloud.point(0).dot(normal) > 0) normal = -normal;
      const Eigen::Vector3d values(1e-3, 1, 1);
      const Eigen::Matrix3d regularized = vectors * values.asDiagonal() * vectors.transpose();
      OrderedNeighbors tree{count};
      small_gicp::estimate_local_features<small_gicp::NormalCovarianceSetter<small_gicp::PointCloud>>(
        cloud, tree, static_cast<int>(count), 0);
      CHECK((cloud.normal(0) - normal).norm() < 1e-11);
      CHECK((cloud.cov(0).topLeftCorner<3, 3>() - regularized).norm() < 1e-11);
      CHECK(cloud.cov(0).row(3).isZero());
      CHECK(cloud.cov(0).col(3).isZero());
    }
  }
}

TEST_CASE("local features retain the insufficient-neighbor and nonfinite fallbacks") {
  std::vector<Eigen::Vector3d> points(5, Eigen::Vector3d(1, 2, 3));
  small_gicp::PointCloud cloud(points);
  OrderedNeighbors tree{4};
  small_gicp::estimate_normals_covariances(cloud, tree, 5);
  for (size_t i = 0; i < cloud.size(); ++i) {
    CHECK(cloud.normal(i).isZero());
    CHECK(cloud.cov(i).topLeftCorner<3, 3>().isIdentity());
    CHECK(cloud.cov(i).row(3).isZero());
    CHECK(cloud.cov(i).col(3).isZero());
  }
  cloud.point(0).x() = std::numeric_limits<double>::quiet_NaN();
  tree.count = 5;
  small_gicp::estimate_local_features<small_gicp::NormalCovarianceSetter<small_gicp::PointCloud>>(
    cloud, tree, 5, 0);
  CHECK(cloud.normal(0).isZero());
  CHECK(cloud.cov(0).topLeftCorner<3, 3>().isIdentity());
}
