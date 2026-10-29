#include <robot_core/test/doctest.h>
#include <small_gicp/ann/incremental_voxelmap.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <array>
#include <type_traits>
#include <vector>

namespace {
template <size_t Rows,size_t Cols>
void checkApprox(const earth::mantle::collection::Matrix<double,Rows,Cols>& actual,
                 const earth::mantle::collection::Matrix<double,Rows,Cols>& expected) {
  for(size_t row=0;row<Rows;++row)
    for(size_t col=0;col<Cols;++col)
      CHECK(actual[row][col]==doctest::Approx(expected[row][col]).epsilon(1e-12));
}
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
  using Point=small_gicp::PointCloud::NativePoint;
  using Covariance=small_gicp::PointCloud::NativeCovariance;
  small_gicp::PointCloud cloud(std::vector<std::array<double,3>>{{1,2,3}});
  cloud.normal(0)=Point::from({0.2,-0.3,0.4,0});
  for(size_t row=0;row<4;++row)
    for(size_t col=0;col<4;++col) cloud.cov(0)[row][col]=10*row+col;
  // Exact 90-degree Z rotation plus translation. Expected values below are
  // evaluated analytically, independently of the production matrix products.
  const auto transform=small_gicp::Transform::from({0,-1,0,0.3, 1,0,0,0.7, 0,0,1,-1, 0,0,0,1});
  const auto point=Point::from({-1.7,1.7,2,1});
  Container container;
  typename Container::Setting setting;
  container.add(setting,point,cloud,0,transform);
  container.finalize();
  REQUIRE(container.size()==1);
  checkApprox(small_gicp::traits::point(container,0),point);
  CHECK(small_gicp::traits::has_normals(container)==HasNormals);
  CHECK(small_gicp::traits::has_covs(container)==HasCovariances);
  if constexpr (HasNormals) {
    checkApprox(small_gicp::traits::normal(container,0),Point::from({0.3,0.2,0.4,0}));
  }
  if constexpr (HasCovariances) {
    const auto covariance=Covariance::from({0.77,-3.17,0.7,-3.1, -14.87,39.27,-1.7,26.1,
                                          7,-17,0,-10, -21.1,53.1,-1,33});
    checkApprox(small_gicp::traits::cov(container,0),covariance);
  }
  const auto snapshot=small_gicp::traits::point(container,0);
  Container copy=container;
  copy.points[0][0][0]+=2;
  const auto original_point=small_gicp::traits::point(container,0);
  for(size_t axis=0;axis<4;++axis) CHECK(original_point[axis][0]==snapshot[axis][0]);
  CHECK(small_gicp::traits::point(copy,0)[0][0]==doctest::Approx(snapshot[0][0]+2));
}
}

TEST_CASE("native flat containers preserve every optional attribute combination") {
  checkContainerAttributes<false,false>();
  checkContainerAttributes<true,false>();
  checkContainerAttributes<false,true>();
  checkContainerAttributes<true,true>();
}

TEST_CASE("native flat container preserves rejection boundaries and neighbor order") {
  using Point=small_gicp::PointCloud::NativePoint;
  small_gicp::FlatContainerNormalCov container;
  small_gicp::FlatContainerNormalCov::Setting setting;
  setting.min_sq_dist_in_cell=0.25;
  setting.max_num_points_in_cell=2;
  small_gicp::PointCloud cloud(std::vector<std::array<double,3>>{{0,0,0}});
  cloud.normal(0)=Point::from({0,0,1,0});
  cloud.cov(0)=earth::mantle::math::identityMatrix<double,4>();
  size_t indices[3]{};
  double distances[3]{};
  CHECK(container.nearest_neighbor_search(cloud.point(0),indices,distances)==0);
  const auto transform=earth::mantle::math::identityMatrix<double,4>();
  container.add(setting,cloud.point(0),cloud,0,transform);
  cloud.point(0)[0][0]=0.25;  // Strictly inside the minimum-distance bound.
  container.add(setting,cloud.point(0),cloud,0,transform);
  CHECK(container.size()==1);
  CHECK(container.normals.size()==1);
  CHECK(container.covs.size()==1);
  cloud.point(0)[0][0]=0.5;  // Equality is accepted by the original strict '<'.
  container.add(setting,cloud.point(0),cloud,0,transform);
  REQUIRE(container.size()==2);
  cloud.point(0)[0][0]=1.5;
  container.add(setting,cloud.point(0),cloud,0,transform);
  CHECK(container.size()==2);
  CHECK(container.normals.size()==2);
  CHECK(container.covs.size()==2);
  const auto query=Point::from({0.4,0,0,1});
  REQUIRE(container.knn_search(query,3,indices,distances)==2);
  CHECK(indices[0]==1);
  CHECK(indices[1]==0);
  CHECK(distances[0]==doctest::Approx(0.01));
  CHECK(distances[1]==doctest::Approx(0.16));
}

TEST_CASE("native flat container remains usable through voxel-map traits") {
  small_gicp::PointCloud cloud(std::vector<std::array<double,3>>{{0.1,0.2,0.3}});
  cloud.cov(0)=earth::mantle::math::identityMatrix<double,4>();
  small_gicp::IncrementalVoxelMap<small_gicp::FlatContainerCov> map(1.0);
  map.insert(cloud);
  cloud.point(0)[0][0]=0.4;
  map.insert(cloud);
  REQUIRE(map.size()==1);
  size_t index=0;
  double distance=0;
  REQUIRE(map.nearest_neighbor_search(cloud.point(0),&index,&distance)==1);
  const auto point=small_gicp::traits::point(map,index);
  const auto covariance=small_gicp::traits::cov(map,index);
  for(size_t row=0;row<4;++row) {
    CHECK(point[row][0]==cloud.point(0)[row][0]);
    for(size_t col=0;col<4;++col) CHECK(covariance[row][col]==cloud.cov(0)[row][col]);
  }
}
