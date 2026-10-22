#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <small_gicp/ros/ros_impl.hpp>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

namespace {
struct Field {
  static constexpr uint8_t FLOAT32=7;
  std::string name;
  uint32_t offset{},count{};
  uint8_t datatype{};
};
struct Message {
  struct Header { std::string frame_id;uint64_t stamp{}; } header;
  uint32_t width{},height{},point_step{},row_step{};
  bool is_bigendian{},is_dense{};
  std::vector<Field> fields;
  std::vector<uint8_t> data;
};
}

TEST_CASE("PointCloud2 conversion keeps XYZ and normal output with native storage") {
  small_gicp::PointCloud input(std::vector<Eigen::Vector3d>{{1,2,3},{-4,5,-6}});
  small_gicp::eigen_view(input.normal(0))=Eigen::Vector4d(0,0,1,0);
  small_gicp::eigen_view(input.normal(1))=Eigen::Vector4d(0,1,0,0);
  const auto message=small_gicp::detail::to_impl<Message,Field,uint64_t>(input,"map",123,true);
  CHECK(message.header.frame_id=="map");
  CHECK(message.header.stamp==123);
  CHECK(message.width==2);
  CHECK(message.fields.size()==6);
  const auto decoded=small_gicp::detail::from_impl<Message,Field>(message);
  REQUIRE(decoded->size()==2);
  CHECK(small_gicp::eigen_view(decoded->point(0))==small_gicp::eigen_view(input.point(0)));
  CHECK(small_gicp::eigen_view(decoded->point(1))==small_gicp::eigen_view(input.point(1)));
  CHECK(decoded->normals.size()==2);
  CHECK(decoded->covs.size()==2);
  float normal_z=0;
  std::memcpy(&normal_z,message.data.data()+message.fields[5].offset,sizeof(normal_z));
  CHECK(normal_z==1);
}

TEST_CASE("PointCloud2 native storage conversion still skips invalid points and accepts empty messages") {
  small_gicp::PointCloud input(std::vector<Eigen::Vector3d>{{1,2,3},{4,5,6}});
  auto message=small_gicp::detail::to_impl<Message,Field,uint64_t>(input,"map",123,false);
  const float invalid=std::numeric_limits<float>::quiet_NaN();
  std::memcpy(message.data.data(),&invalid,sizeof(invalid));
  const auto decoded=small_gicp::detail::from_impl<Message,Field>(message);
  REQUIRE(decoded->size()==1);
  CHECK(small_gicp::eigen_view(decoded->point(0))==small_gicp::eigen_view(input.point(1)));
  CHECK(small_gicp::detail::from_impl<Message,Field>(Message{})->empty());
}
