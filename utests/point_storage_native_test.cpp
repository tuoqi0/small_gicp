#include <robot_core/test/doctest.h>
#include <small_gicp/points/point_cloud.hpp>
#include <math/fixed_matrix_ops.h>
#include <array>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
using Cloud=small_gicp::PointCloud;
static_assert(std::is_same<decltype(std::declval<Cloud&>().point(0)), Cloud::NativePoint&>::value,
    "production point access must return native storage");
static_assert(std::is_same<decltype(std::declval<const Cloud&>().cov(0)), const Cloud::NativeCovariance&>::value,
    "const covariance access must remain a native reference");
static_assert(std::is_same<decltype(Cloud::points)::value_type,
    earth::mantle::collection::Matrix<double,4,1>>::value,"point storage must be native");
static_assert(std::is_same<decltype(Cloud::normals)::value_type,
    earth::mantle::collection::Matrix<double,4,1>>::value,"normal storage must be native");
static_assert(std::is_same<decltype(Cloud::covs)::value_type,
    earth::mantle::collection::Matrix<double,4,4>>::value,"covariance storage must be native");
template <size_t Rows,size_t Cols>
void checkSame(const earth::mantle::collection::Matrix<double,Rows,Cols>& actual,
               const earth::mantle::collection::Matrix<double,Rows,Cols>& expected) {
  for(size_t row=0;row<Rows;++row)
    for(size_t col=0;col<Cols;++col) CHECK(actual[row][col]==expected[row][col]);
}
}

TEST_CASE("native point buffers preserve constructor values and mutable references") {
  Cloud cloud(std::vector<std::array<double,3>>{{1,2,3},{4,5,6}});
  REQUIRE(cloud.size()==2);
  checkSame(cloud.point(0),Cloud::NativePoint::from({1,2,3,1}));
  auto& point=cloud.point(0);
  point[1][0]=8;
  CHECK(cloud.points[0][1][0]==8);
  cloud.normal(0)=Cloud::NativePoint::from({0,0,1,0});
  CHECK(cloud.normals[0][2][0]==1);
  const auto snapshot=small_gicp::traits::point(cloud,0);
  const Cloud& read_only=cloud;
  const auto& read_point=read_only.point(0);
  cloud.point(0)[0][0]=10;
  CHECK(snapshot[0][0]==1);
  CHECK(cloud.point(0)[0][0]==10);
  CHECK(read_point[0][0]==10);
}

TEST_CASE("native covariance references preserve row-major storage and const value semantics") {
  Cloud cloud;cloud.resize(1);
  Cloud::NativeCovariance covariance{};
  for(size_t row=0;row<4;++row)
    for(size_t col=0;col<4;++col) covariance[row][col]=10*row+col;
  cloud.cov(0)=covariance;
  CHECK(cloud.covs[0][0][1]==1);
  CHECK(cloud.covs[0][1][0]==10);
  const Cloud& read_only=cloud;
  checkSame(read_only.cov(0),covariance);
  const auto snapshot=small_gicp::traits::cov(cloud,0);
  const auto& read_cov=read_only.cov(0);
  cloud.cov(0)[0][1]=99;
  CHECK(snapshot[0][1]==1);
  CHECK(cloud.covs[0][0][1]==99);
  CHECK(read_cov[0][1]==99);
  const auto before=cloud.cov(0);
  small_gicp::traits::set_cov(cloud,0,earth::mantle::math::transpose(cloud.cov(0)));
  checkSame(cloud.cov(0),earth::mantle::math::transpose(before));
}

TEST_CASE("native point-cloud copies and moves own independent buffers") {
  Cloud original(std::vector<std::array<double,3>>{{1,2,3}});
  original.normal(0)=Cloud::NativePoint::from({0,0,1,0});
  original.cov(0)=earth::mantle::math::identityMatrix<double,4>();
  Cloud copy=original;
  copy.point(0)[0][0]=20;
  copy.normal(0)[2][0]=0;
  copy.cov(0)[1][1]=7;
  CHECK(original.point(0)[0][0]==1);
  CHECK(original.normal(0)[2][0]==1);
  CHECK(original.cov(0)[1][1]==1);
  const auto* point_buffer=copy.points.data();
  const auto* normal_buffer=copy.normals.data();
  const auto* covariance_buffer=copy.covs.data();
  Cloud moved=std::move(copy);
  CHECK(moved.points.data()==point_buffer);
  CHECK(moved.normals.data()==normal_buffer);
  CHECK(moved.covs.data()==covariance_buffer);
  CHECK(moved.point(0)[0][0]==20);
  CHECK(moved.cov(0)[1][1]==7);
  moved.resize(2);
  CHECK(moved.point(0)[0][0]==20);
  checkSame(moved.point(1),Cloud::NativePoint{});
  checkSame(moved.normal(1),Cloud::NativePoint{});
  checkSame(moved.cov(1),Cloud::NativeCovariance{});
  Cloud assigned;
  const auto* moved_buffer=moved.points.data();
  assigned=std::move(moved);
  CHECK(assigned.points.data()==moved_buffer);
  CHECK(assigned.point(0)[0][0]==20);
}
