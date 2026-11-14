#include <random>
#include <small_gicp/util/fast_floor.hpp>
#include <small_gicp/util/vector3i_hash.hpp>

#include <robot_core/test/doctest.h>
#include <cmath>

using namespace small_gicp;

TEST_CASE("VectorTest.HashTest") {
  std::mt19937 mt;

  for (int i = 0; i < 1000; i++) {
    std::uniform_int_distribution<> dist(-1000, 1000);
    const VoxelCoord v{{dist(mt), dist(mt), dist(mt)}};
    CHECK((XORVector3iHash::hash(v)) == (XORVector3iHash()(v)));
    CHECK((XORVector3iHash::equal(v, v)));
  }
}

TEST_CASE("VectorTest.FloorTest") {
  std::mt19937 mt;
  for (int i = 0; i < 1000; i++) {
    std::uniform_real_distribution<> dist(-1000.0, 1000.0);
    const std::array<double, 4> v{{dist(mt), dist(mt), dist(mt), 1.0}};
    const auto floor1 = fast_floor(v);
    CHECK((floor1[0]) == (std::floor(v[0])));
    CHECK((floor1[1]) == (std::floor(v[1])));
    CHECK((floor1[2]) == (std::floor(v[2])));
    CHECK((floor1[3]) == (std::floor(v[3])));
  }
}
