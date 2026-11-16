#include <robot_core/test/doctest.h>
#include <small_gicp/ros/ros_impl.hpp>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

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
void checkPoint(const small_gicp::PointCloud::NativePoint& actual,
                const small_gicp::PointCloud::NativePoint& expected) {
  for(size_t axis=0;axis<4;++axis) CHECK(actual[axis][0]==expected[axis][0]);
}
}

TEST_CASE("PointCloud2 conversion keeps XYZ and normal output with native storage") {
  small_gicp::PointCloud input(std::vector<std::array<double,3>>{{1,2,3},{-4,5,-6}});
  input.normal(0)=small_gicp::PointCloud::NativePoint::from({0,0,1,0});
  input.normal(1)=small_gicp::PointCloud::NativePoint::from({0,1,0,0});
  const auto message=small_gicp::detail::to_impl<Message,Field,uint64_t>(input,"map",123,true);
  CHECK(message.header.frame_id=="map");
  CHECK(message.header.stamp==123);
  CHECK(message.width==2);
  CHECK(message.fields.size()==6);
  const auto decoded=small_gicp::detail::from_impl<Message,Field>(message);
  REQUIRE(decoded->size()==2);
  checkPoint(decoded->point(0),input.point(0));
  checkPoint(decoded->point(1),input.point(1));
  CHECK(decoded->normals.size()==2);
  CHECK(decoded->covs.size()==2);
  float normal_z=0;
  std::memcpy(&normal_z,message.data.data()+message.fields[5].offset,sizeof(normal_z));
  CHECK(normal_z==1);
}

TEST_CASE("PointCloud2 native storage conversion still skips invalid points and accepts empty messages") {
  small_gicp::PointCloud input(std::vector<std::array<double,3>>{{1,2,3},{4,5,6}});
  auto message=small_gicp::detail::to_impl<Message,Field,uint64_t>(input,"map",123,false);
  const float invalid=std::numeric_limits<float>::quiet_NaN();
  std::memcpy(message.data.data(),&invalid,sizeof(invalid));
  const auto decoded=small_gicp::detail::from_impl<Message,Field>(message);
  REQUIRE(decoded->size()==1);
  checkPoint(decoded->point(0),input.point(1));
  CHECK(small_gicp::detail::from_impl<Message,Field>(Message{})->empty());
}

namespace {
void put_le(std::vector<uint8_t>& data, size_t at, float value) {
  std::memcpy(data.data()+at, &value, sizeof(value));
}
void put_be(std::vector<uint8_t>& data, size_t at, float value) {
  uint8_t b[sizeof(value)]; std::memcpy(b, &value, sizeof(value));
  for (size_t i=0;i<sizeof(value);++i) data[at+i]=b[sizeof(value)-1-i];
}
Message organized() {
  Message m; m.width=2; m.height=2; m.point_step=16; m.row_step=40; m.data.assign(80,0xa5);
  m.fields={{"x",0,1,Field::FLOAT32},{"y",4,1,Field::FLOAT32},{"z",8,1,Field::FLOAT32}};
  for(size_t row=0;row<2;++row) for(size_t col=0;col<2;++col) {
    size_t at=row*m.row_step+col*m.point_step; put_le(m.data,at+0,float(1+row*2+col)); put_le(m.data,at+4,float(2)); put_le(m.data,at+8,float(3));
  }
  return m;
}
}
TEST_CASE("PointCloud2 parser honors organized row padding") {
  auto m=organized(); auto decoded=small_gicp::detail::from_impl<Message,Field>(m);
  REQUIRE(decoded->size()==4); CHECK(decoded->point(2)[0][0]==3); CHECK(decoded->point(3)[0][0]==4);
}
TEST_CASE("PointCloud2 parser decodes big endian and rejects malformed layouts") {
  Message m; m.width=1; m.height=1; m.point_step=12; m.row_step=12; m.is_bigendian=true; m.data.assign(12,0);
  m.fields={{"x",0,1,Field::FLOAT32},{"y",4,1,Field::FLOAT32},{"z",8,1,Field::FLOAT32}};
  put_be(m.data,0,1.25f); put_be(m.data,4,-2.5f); put_be(m.data,8,3.75f);
  auto decoded=small_gicp::detail::from_impl<Message,Field>(m); REQUIRE(decoded->size()==1);
  CHECK(decoded->point(0)[0][0]==1.25); CHECK(decoded->point(0)[1][0]==-2.5); CHECK(decoded->point(0)[2][0]==3.75);
  for(int damage=0; damage<4; ++damage) { auto bad=m; if(damage==0) bad.fields[0].datatype=1; if(damage==1) bad.fields[0].count=2; if(damage==2) bad.row_step=8; if(damage==3) bad.data.resize(11); CHECK(small_gicp::detail::from_impl<Message,Field>(bad)->empty()); }
}
