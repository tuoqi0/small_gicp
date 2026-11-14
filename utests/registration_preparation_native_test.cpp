#include <robot_core/test/doctest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <math/fixed_matrix_ops.h>
#include <small_gicp/ann/kdtree_omp.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/registration/reduction_omp.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/registration/registration_preparation.hpp>
#include <vector>

namespace {
using Cloud = small_gicp::PointCloud;

std::uint64_t bits(double value)
{
  std::uint64_t result;
  static_assert(sizeof(result) == sizeof(value), "double bit contract requires 64-bit storage");
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::vector<std::array<float, 3>> wall()
{
  std::vector<std::array<float, 3>> points;
  for (int y = -4; y <= 4; ++y)
    for (int z = -4; z <= 4; ++z) points.push_back({{4, y * .2F, z * .2F}});
  return points;
}
}  // namespace

TEST_CASE("native registration preparation search and solve need no PCL or Eigen")
{
  REQUIRE(small_gicp::prepare_registration_cloud({}, .25, 20, 2, 1.0)->empty());
  std::vector<std::array<float, 3>> points;
  for (int x = -3; x <= 3; ++x)
    for (int y = -3; y <= 3; ++y)
      for (int z = -3; z <= 3; ++z) points.push_back({{x * .5F, y * .5F, z * .5F}});
  const auto cloud = small_gicp::prepare_registration_cloud(points, .25, 20, 2, 1.0);
  REQUIRE(cloud->size() == points.size());
  Cloud::Ptr nonnull_empty = small_gicp::prepare_registration_cloud({}, .25, 20, 2, 1.0);
  REQUIRE(nonnull_empty != nullptr);
  CHECK(nonnull_empty->points.empty());
  CHECK(nonnull_empty->normals.empty());
  CHECK(nonnull_empty->covs.empty());
  small_gicp::KdTree<Cloud> tree(cloud, small_gicp::KdTreeBuilderOMP(2));
  for (std::size_t i = 0; i < cloud->size(); ++i) {
    std::size_t index = 0;
    double squared_distance = 1;
    REQUIRE(tree.nearest_neighbor_search(cloud->point(i), &index, &squared_distance));
    CHECK(index == i);
    CHECK(squared_distance == 0);
    CHECK(cloud->point(i)[3][0] == 1);
    for (std::size_t j = 0; j < 16; ++j) CHECK(std::isfinite(cloud->cov(i).raw_pointer()[j]));
  }
  small_gicp::Registration<small_gicp::GICPFactor, small_gicp::ParallelReductionOMP> registration;
  registration.reduction.num_threads = 1;
  const auto identity = earth::mantle::math::identityMatrix<double, 4>();
  const auto result = registration.align(*cloud, *cloud, tree, identity);
  REQUIRE(result.converged);
  CHECK(result.num_inliers == cloud->size());
  CHECK(result.error == 0);
  for (std::size_t i = 0; i < 16; ++i)
    CHECK(result.T_target_source.raw_pointer()[i] == doctest::Approx(identity.raw_pointer()[i]));
}

TEST_CASE("geometric normals retain their physical radius separately from GICP covariance")
{
  std::vector<std::array<float, 3>> points;
  for (int y = -10; y <= 10; ++y)
    for (int z = -10; z <= 10; ++z) points.push_back({{0, y * .2F, z * .2F}});
  const auto cloud = small_gicp::prepare_registration_cloud(points, .25, 10, 2, 1.0);
  REQUIRE_FALSE(cloud->empty());
  REQUIRE(cloud->normals.size() == cloud->size());
  for (std::size_t i = 0; i < cloud->size(); ++i) {
    CHECK(std::abs(cloud->normal(i)[0][0]) == doctest::Approx(1));
    CHECK(cloud->normal(i)[1][0] == doctest::Approx(0));
    CHECK(cloud->normal(i)[2][0] == doctest::Approx(0));
    CHECK(cloud->normal(i)[3][0] == 0);
    for (std::size_t j = 0; j < 16; ++j) CHECK(std::isfinite(cloud->cov(i).raw_pointer()[j]));
  }
}

TEST_CASE("voxel centroids return to float precision without changing the input")
{
  const std::vector<std::array<float, 3>> points{{{.1F, 2.25F, 3.25F}}, {{.2F, 2.25F, 3.25F}}};
  const auto before = points;
  const auto cloud = small_gicp::prepare_registration_cloud(points, 1.0, 20, 1, .5);
  REQUIRE(cloud != nullptr);
  REQUIRE(cloud->size() == 1);
  CHECK(points == before);
  CHECK(bits(cloud->point(0)[0][0]) == bits(static_cast<double>(.15F)));
  CHECK(cloud->point(0)[1][0] == 2.25);
  CHECK(cloud->point(0)[2][0] == 3.25);
  CHECK(cloud->point(0)[3][0] == 1);
  const double unrounded_centroid = (static_cast<double>(.1F) + static_cast<double>(.2F)) / 2;
  CHECK(bits(cloud->point(0)[0][0]) != bits(unrounded_centroid));
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t column = 0; column < 4; ++column)
      CHECK(cloud->cov(0)[row][column] == (row == column && row < 3 ? 1.0 : 0.0));
  for (std::size_t axis = 0; axis < 3; ++axis) CHECK(std::isnan(cloud->normal(0)[axis][0]));
  CHECK(cloud->normal(0)[3][0] == 0);
}

TEST_CASE("oblique geometry preserves float covariance storage and origin facing normals")
{
  std::vector<std::array<float, 3>> points;
  for (int x = -2; x <= 2; ++x)
    for (int y = -2; y <= 2; ++y)
      points.push_back({{static_cast<float>(x), static_cast<float>(y), static_cast<float>(2 - x - y)}});
  const auto before = points;
  const auto cloud = small_gicp::prepare_registration_cloud(points, .1, 25, 2, 20.0);
  REQUIRE(cloud->size() == points.size());
  CHECK(points == before);
  for (std::size_t i = 0; i < cloud->size(); ++i) {
    for (std::size_t row = 0; row < 4; ++row) {
      for (std::size_t column = 0; column < 4; ++column) {
        const double value = cloud->cov(i)[row][column];
        REQUIRE(std::isfinite(value));
        CHECK(bits(value) == bits(static_cast<double>(static_cast<float>(value))));
        if (row == 3 || column == 3) CHECK(value == 0);
        else CHECK(value == doctest::Approx(row == column ? .667 : -.333).epsilon(1e-6));
      }
    }
    for (std::size_t axis = 0; axis < 3; ++axis) {
      const double normal = cloud->normal(i)[axis][0];
      REQUIRE(std::isfinite(normal));
      CHECK(normal == doctest::Approx(-.5773502691896258).epsilon(1e-6));
      CHECK(bits(normal) == bits(static_cast<double>(static_cast<float>(normal))));
    }
    CHECK(cloud->normal(i)[3][0] == 0);
  }
}

TEST_CASE("radius support changes geometric normals without changing GICP weights")
{
  const auto points = wall();
  const auto narrow = small_gicp::prepare_registration_cloud(points, .05, 20, 2, .1);
  const auto supported = small_gicp::prepare_registration_cloud(points, .05, 20, 2, .5);
  const auto sparse_knn = small_gicp::prepare_registration_cloud(points, .05, 4, 2, .5);
  REQUIRE(narrow->size() == points.size());
  REQUIRE(supported->size() == narrow->size());
  REQUIRE(sparse_knn->size() == narrow->size());
  for (std::size_t i = 0; i < narrow->size(); ++i) {
    for (std::size_t axis = 0; axis < 4; ++axis) {
      CHECK(bits(narrow->point(i)[axis][0]) == bits(supported->point(i)[axis][0]));
      CHECK(bits(sparse_knn->normal(i)[axis][0]) == bits(supported->normal(i)[axis][0]));
    }
    for (std::size_t j = 0; j < 16; ++j)
      CHECK(bits(narrow->cov(i).raw_pointer()[j]) == bits(supported->cov(i).raw_pointer()[j]));
    for (std::size_t axis = 0; axis < 3; ++axis) CHECK(std::isnan(narrow->normal(i)[axis][0]));
    CHECK(supported->normal(i)[0][0] == doctest::Approx(-1));
    CHECK(supported->normal(i)[1][0] == 0);
    CHECK(supported->normal(i)[2][0] == 0);
    CHECK(supported->normal(i)[3][0] == 0);
    CHECK(supported->cov(i)[0][0] == doctest::Approx(.001).epsilon(1e-6));
    CHECK(sparse_knn->cov(i)[0][0] == 1);
  }
}

TEST_CASE("collinear radius neighborhoods keep invalid normals and finite covariance")
{
  std::vector<std::array<float, 3>> points;
  for (int i = -3; i <= 3; ++i) points.push_back({{2, i * .3F, 0}});
  const auto cloud = small_gicp::prepare_registration_cloud(points, .05, 7, 2, 3.0);
  REQUIRE(cloud->size() == points.size());
  for (std::size_t i = 0; i < cloud->size(); ++i) {
    for (std::size_t axis = 0; axis < 3; ++axis) CHECK(std::isnan(cloud->normal(i)[axis][0]));
    CHECK(cloud->normal(i)[3][0] == 0);
    for (std::size_t j = 0; j < 16; ++j) CHECK(std::isfinite(cloud->cov(i).raw_pointer()[j]));
  }
}

TEST_CASE("invalid radius retains the prepared cloud and untouched normal buffers")
{
  const auto points = wall();
  for (const double radius : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::infinity()}) {
    const auto cloud = small_gicp::prepare_registration_cloud(points, .05, 20, 1, radius);
    REQUIRE(cloud != nullptr);
    REQUIRE(cloud->size() == points.size());
    REQUIRE(cloud->normals.size() == cloud->size());
    REQUIRE(cloud->covs.size() == cloud->size());
    for (std::size_t i = 0; i < cloud->size(); ++i) {
      CHECK(cloud->point(i)[3][0] == 1);
      for (std::size_t axis = 0; axis < 4; ++axis) CHECK(cloud->normal(i)[axis][0] == 0);
      for (std::size_t j = 0; j < 16; ++j) CHECK(std::isfinite(cloud->cov(i).raw_pointer()[j]));
    }
  }
}
