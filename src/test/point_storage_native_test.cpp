#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <small_gicp/points/point_cloud.hpp>
#include <type_traits>
#include <utility>

namespace {
using Cloud=small_gicp::PointCloud;
static_assert(std::is_same<decltype(std::declval<Cloud&>().point(0)), Cloud::NativePoint&>::value,
    "production point access must return native storage, not an Eigen view");
static_assert(std::is_same<decltype(std::declval<const Cloud&>().cov(0)), const Cloud::NativeCovariance&>::value,
    "const covariance access must remain a native reference");
static_assert(std::is_same<decltype(Cloud::points)::value_type,
    earth::mantle::collection::Matrix<double,4,1>>::value,"point storage must be native");
static_assert(std::is_same<decltype(Cloud::normals)::value_type,
    earth::mantle::collection::Matrix<double,4,1>>::value,"normal storage must be native");
static_assert(std::is_same<decltype(Cloud::covs)::value_type,
    earth::mantle::collection::Matrix<double,4,4>>::value,"covariance storage must be native");
}

TEST_CASE("native point buffers preserve constructor values and mutable views") {
  Cloud cloud(std::vector<Eigen::Vector3d>{{1,2,3},{4,5,6}});
  REQUIRE(cloud.size()==2);
  CHECK(small_gicp::eigen_view(cloud.point(0))==Eigen::Vector4d(1,2,3,1));
  auto point=small_gicp::eigen_view(cloud.point(0));
  point[1]=8;
  CHECK(cloud.points[0][1][0]==8);
  small_gicp::eigen_view(cloud.normal(0))=Eigen::Vector4d(0,0,1,0);
  CHECK(cloud.normals[0][2][0]==1);
  const auto snapshot=small_gicp::eigen_matrix(small_gicp::traits::point(cloud,0));
  const Cloud& read_only=cloud;
  auto head=small_gicp::eigen_view(read_only.point(0)).head<3>();
  small_gicp::eigen_view(cloud.point(0))[0]=10;
  CHECK(snapshot[0]==1);
  CHECK(small_gicp::eigen_view(cloud.point(0))[0]==10);
  CHECK(head[0]==10);
}

TEST_CASE("native covariance views preserve row-major storage and const value semantics") {
  Cloud cloud;cloud.resize(1);
  Eigen::Matrix4d covariance;
  for(int row=0;row<4;++row)
    for(int col=0;col<4;++col) covariance(row,col)=10*row+col;
  small_gicp::eigen_view(cloud.cov(0))=covariance;
  CHECK(cloud.covs[0][0][1]==1);
  CHECK(cloud.covs[0][1][0]==10);
  const Cloud& read_only=cloud;
  CHECK(small_gicp::eigen_view(read_only.cov(0))==covariance);
  const auto snapshot=small_gicp::eigen_matrix(small_gicp::traits::cov(cloud,0));
  auto corner=small_gicp::eigen_view(read_only.cov(0)).topLeftCorner<2,2>();
  small_gicp::eigen_view(cloud.cov(0))(0,1)=99;
  CHECK(snapshot(0,1)==1);
  CHECK(cloud.covs[0][0][1]==99);
  CHECK(corner(0,1)==99);
  const Eigen::Matrix4d before = small_gicp::eigen_view(cloud.cov(0));
  small_gicp::traits::set_cov(cloud, 0, small_gicp::eigen_view(cloud.cov(0)).transpose());
  CHECK(small_gicp::eigen_view(cloud.cov(0)) == before.transpose());
}

TEST_CASE("native point-cloud copies and moves own independent buffers") {
  Cloud original(std::vector<Eigen::Vector3d>{{1,2,3}});
  small_gicp::eigen_view(original.normal(0))=Eigen::Vector4d(0,0,1,0);
  small_gicp::eigen_view(original.cov(0))=Eigen::Matrix4d::Identity();
  Cloud copy=original;
  small_gicp::eigen_view(copy.point(0))[0]=20;
  small_gicp::eigen_view(copy.normal(0))[2]=0;
  small_gicp::eigen_view(copy.cov(0))(1,1)=7;
  CHECK(small_gicp::eigen_view(original.point(0))[0]==1);
  CHECK(small_gicp::eigen_view(original.normal(0))[2]==1);
  CHECK(small_gicp::eigen_view(original.cov(0))(1,1)==1);
  const auto* point_buffer=copy.points.data();
  const auto* normal_buffer=copy.normals.data();
  const auto* covariance_buffer=copy.covs.data();
  Cloud moved=std::move(copy);
  CHECK(moved.points.data()==point_buffer);
  CHECK(moved.normals.data()==normal_buffer);
  CHECK(moved.covs.data()==covariance_buffer);
  CHECK(small_gicp::eigen_view(moved.point(0))[0]==20);
  CHECK(small_gicp::eigen_view(moved.cov(0))(1,1)==7);
  moved.resize(2);
  CHECK(small_gicp::eigen_view(moved.point(0))[0]==20);
  CHECK(small_gicp::eigen_view(moved.point(1)).isZero());
  CHECK(small_gicp::eigen_view(moved.normal(1)).isZero());
  CHECK(small_gicp::eigen_view(moved.cov(1)).isZero());
  Cloud assigned;
  const auto* moved_buffer=moved.points.data();
  assigned=std::move(moved);
  CHECK(assigned.points.data()==moved_buffer);
  CHECK(small_gicp::eigen_view(assigned.point(0))[0]==20);
}
