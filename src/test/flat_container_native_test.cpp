#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <small_gicp/ann/incremental_voxelmap.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <type_traits>

namespace {
template <bool HasNormals, bool HasCovariances>
void checkContainerAttributes() {
  using Container = small_gicp::FlatContainer<HasNormals, HasCovariances>;
  static_assert(std::is_same<typename decltype(Container::points)::value_type,
      earth::mantle::collection::Matrix<double, 4, 1>>::value,
      "flat points must own native storage");
  if constexpr (HasNormals) {
    static_assert(std::is_same<typename decltype(Container::normals)::value_type,
        earth::mantle::collection::Matrix<double, 4, 1>>::value,
        "flat normals must own native storage");
  }
  if constexpr (HasCovariances) {
    static_assert(std::is_same<typename decltype(Container::covs)::value_type,
        earth::mantle::collection::Matrix<double, 4, 4>>::value,
        "flat covariances must own native storage");
  }
  small_gicp::PointCloud cloud(std::vector<Eigen::Vector3d>{{1, 2, 3}});
  small_gicp::eigen_view(cloud.normal(0)) = Eigen::Vector4d(0.2, -0.3, 0.4, 0);
  for (int row = 0; row < 4; ++row)
    for (int col = 0; col < 4; ++col) small_gicp::eigen_view(cloud.cov(0))(row, col) = 10 * row + col;
  Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
  transform.translation() = Eigen::Vector3d(0.3, 0.7, -1);
  transform.linear() = Eigen::AngleAxisd(0.4,
      Eigen::Vector3d(1, -2, 3).normalized()).toRotationMatrix();
  Container container;
  typename Container::Setting setting;
  const Eigen::Vector4d point = transform * small_gicp::eigen_view(cloud.point(0));
  container.add(setting, small_gicp::native_matrix(point), cloud, 0, small_gicp::native_transform(transform));
  container.finalize();
  REQUIRE(container.size() == 1);
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(container, 0)).isApprox(point, 1e-12));
  CHECK(small_gicp::traits::has_normals(container) == HasNormals);
  CHECK(small_gicp::traits::has_covs(container) == HasCovariances);
  if constexpr (HasNormals) {
    const Eigen::Vector4d normal = transform.matrix() * small_gicp::eigen_view(cloud.normal(0));
    CHECK(small_gicp::eigen_matrix(small_gicp::traits::normal(container, 0)).isApprox(normal, 1e-12));
  }
  if constexpr (HasCovariances) {
    const Eigen::Matrix4d covariance =
        transform.matrix() * small_gicp::eigen_view(cloud.cov(0)) * transform.matrix().transpose();
    CHECK(small_gicp::eigen_matrix(small_gicp::traits::cov(container, 0)).isApprox(covariance, 1e-12));
  }
  const auto snapshot = small_gicp::eigen_matrix(small_gicp::traits::point(container, 0));
  Container copy = container;
  copy.points[0][0][0] += 2;
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(container, 0)) == snapshot);
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(copy, 0))[0] == doctest::Approx(snapshot[0] + 2));
}
}

TEST_CASE("native flat containers preserve every optional attribute combination") {
  checkContainerAttributes<false, false>();
  checkContainerAttributes<true, false>();
  checkContainerAttributes<false, true>();
  checkContainerAttributes<true, true>();
}

TEST_CASE("native flat container preserves rejection boundaries and neighbor order") {
  small_gicp::FlatContainerNormalCov container;
  small_gicp::FlatContainerNormalCov::Setting setting;
  setting.min_sq_dist_in_cell = 0.25;
  setting.max_num_points_in_cell = 2;
  small_gicp::PointCloud cloud(std::vector<Eigen::Vector3d>{{0, 0, 0}});
  small_gicp::eigen_view(cloud.normal(0)) = Eigen::Vector4d(0, 0, 1, 0);
  small_gicp::eigen_view(cloud.cov(0)).setIdentity();
  size_t indices[3]{};
  double distances[3]{};
  CHECK(container.nearest_neighbor_search(small_gicp::eigen_view(cloud.point(0)), indices, distances) == 0);
  const Eigen::Isometry3d transform = Eigen::Isometry3d::Identity();
  container.add(setting, cloud.point(0), cloud, 0, small_gicp::native_transform(transform));
  small_gicp::eigen_view(cloud.point(0))[0] = 0.25;  // Strictly inside the minimum-distance bound.
  container.add(setting, cloud.point(0), cloud, 0, small_gicp::native_transform(transform));
  CHECK(container.size() == 1);
  CHECK(container.normals.size() == 1);
  CHECK(container.covs.size() == 1);
  small_gicp::eigen_view(cloud.point(0))[0] = 0.5;  // Equality is accepted by the original strict '<'.
  container.add(setting, cloud.point(0), cloud, 0, small_gicp::native_transform(transform));
  REQUIRE(container.size() == 2);
  small_gicp::eigen_view(cloud.point(0))[0] = 1.5;
  container.add(setting, cloud.point(0), cloud, 0, small_gicp::native_transform(transform));
  CHECK(container.size() == 2);
  CHECK(container.normals.size() == 2);
  CHECK(container.covs.size() == 2);
  const Eigen::Vector4d query(0.4, 0, 0, 1);
  REQUIRE(container.knn_search(query, 3, indices, distances) == 2);
  CHECK(indices[0] == 1);
  CHECK(indices[1] == 0);
  CHECK(distances[0] == doctest::Approx(0.01));
  CHECK(distances[1] == doctest::Approx(0.16));
}

TEST_CASE("native flat container remains usable through voxel-map traits") {
  small_gicp::PointCloud cloud(std::vector<Eigen::Vector3d>{{0.1, 0.2, 0.3}});
  small_gicp::eigen_view(cloud.cov(0)).setIdentity();
  small_gicp::IncrementalVoxelMap<small_gicp::FlatContainerCov> map(1.0);
  map.insert(cloud);
  small_gicp::eigen_view(cloud.point(0))[0] = 0.4;
  map.insert(cloud);
  REQUIRE(map.size() == 1);
  size_t index = 0;
  double distance = 0;
  REQUIRE(map.nearest_neighbor_search(small_gicp::eigen_view(cloud.point(0)), &index, &distance) == 1);
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::point(map, index)) == small_gicp::eigen_view(cloud.point(0)));
  CHECK(small_gicp::eigen_matrix(small_gicp::traits::cov(map, index)) == small_gicp::eigen_view(cloud.cov(0)));
}
