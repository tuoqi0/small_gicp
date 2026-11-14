#include <robot_core/test/doctest.h>
#include <small_gicp/registration/loop_registration.hpp>
#include <math/fixed_matrix_ops.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <omp.h>

namespace {
namespace alg = earth::mantle::alg;
namespace math = earth::mantle::math;
using Cloud = alg::PointCloudXYZI;
using Transform = alg::loop_closure::RegistrationTransform;

constexpr float kVoxelMetres = .2f;
constexpr int kCovarianceNeighbors = 20;
constexpr int kWorkerThreads = 1;
constexpr float kMaxCorrespondenceMetres = 10.f;
constexpr int kIterations = 100;
constexpr float kHalfSourceRatio = .5f;

Cloud::Ptr lattice() {
  auto cloud = std::make_shared<Cloud>();
  for (int x = -2; x <= 2; ++x) {
    for (int y = -2; y <= 2; ++y) {
      for (int z = -2; z <= 2; ++z) {
        alg::PointXYZINormal point{};
        point.x = static_cast<float>(x);
        point.y = static_cast<float>(y);
        point.z = static_cast<float>(z);
        point.intensity = static_cast<float>(cloud->size());
        cloud->push_back(point);
      }
    }
  }
  return cloud;
}

Cloud::Ptr halfSupportedSource(const Cloud& target) {
  auto source = std::make_shared<Cloud>(target);
  for (auto point : target.points) {
    // Both components have unique .2m voxels. Only the first125 points are
    // within10m of the target, so the sampled source denominator stays250.
    point.x += 64.f;
    source->push_back(point);
  }
  return source;
}

Transform outputSentinel() {
  Transform result;
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 4; ++col)
      result[row][col] = 10.f * static_cast<float>(row) + static_cast<float>(col) + .25f;
  return result;
}

void checkTransformBits(const Transform& actual, const Transform& expected) {
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      CAPTURE(row);
      CAPTURE(col);
      CHECK(std::memcmp(&actual[row][col], &expected[row][col], sizeof(float)) == 0);
    }
  }
}

std::array<float, 8> tuple(const alg::PointXYZINormal& point) {
  return {{point.x, point.y, point.z, point.intensity,
           point.normal_x, point.normal_y, point.normal_z, point.curvature}};
}

void checkCloudUnchanged(const Cloud& actual, const Cloud& expected) {
  REQUIRE(actual.size() == expected.size());
  CHECK(actual.width == expected.width);
  CHECK(actual.height == expected.height);
  CHECK(actual.is_dense == expected.is_dense);
  for (std::size_t i = 0; i < actual.size(); ++i) {
    CAPTURE(i);
    const auto a = tuple(actual[i]);
    const auto b = tuple(expected[i]);
    CHECK(std::memcmp(a.data(), b.data(), sizeof(float) * a.size()) == 0);
  }
}
}  // namespace

TEST_CASE("actual loop GICP provider accepts exactly half of the sampled source") {
  const auto target = lattice();
  const auto source = halfSupportedSource(*target);
  const auto original_source = *source;
  const auto original_target = *target;
  REQUIRE(target->size() == 125);
  REQUIRE(source->size() == 250);
  const auto identity = math::identityMatrix<float, 4>();
  auto output = outputSentinel();
  float ratio = -.25f;
  const int original_workers = omp_get_max_threads();

  const bool accepted = small_gicp::register_loop_gicp(
    source, target, identity, kVoxelMetres, kCovarianceNeighbors, kWorkerThreads,
    kMaxCorrespondenceMetres, kIterations, kHalfSourceRatio, output, ratio);

  REQUIRE(accepted);
  CHECK(ratio == kHalfSourceRatio);
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 4; ++col)
      CHECK(output[row][col] == doctest::Approx(identity[row][col]).epsilon(1e-5));
  CHECK(omp_get_max_threads() == original_workers);
  checkCloudUnchanged(*source, original_source);
  checkCloudUnchanged(*target, original_target);
}

TEST_CASE("actual loop GICP provider rejects insufficient ratio without writing a transform") {
  const auto target = lattice();
  const auto source = halfSupportedSource(*target);
  const auto original_source = *source;
  const auto original_target = *target;
  const auto identity = math::identityMatrix<float, 4>();
  const auto sentinel = outputSentinel();
  auto output = sentinel;
  float ratio = -.25f;
  const float required_ratio = std::nextafter(
    kHalfSourceRatio, std::numeric_limits<float>::infinity());
  const int original_workers = omp_get_max_threads();

  const bool accepted = small_gicp::register_loop_gicp(
    source, target, identity, kVoxelMetres, kCovarianceNeighbors, kWorkerThreads,
    kMaxCorrespondenceMetres, kIterations, required_ratio, output, ratio);

  CHECK_FALSE(accepted);
  CHECK(ratio == kHalfSourceRatio);
  checkTransformBits(output, sentinel);
  CHECK(omp_get_max_threads() == original_workers);
  checkCloudUnchanged(*source, original_source);
  checkCloudUnchanged(*target, original_target);
}

TEST_CASE("actual loop GICP provider keeps both outputs when optimization does not converge") {
  const auto source = lattice();
  const auto target = std::make_shared<Cloud>(*source);
  const auto original_source = *source;
  const auto original_target = *target;
  const auto identity = math::identityMatrix<float, 4>();
  const auto sentinel = outputSentinel();
  auto output = sentinel;
  constexpr float kRatioSentinel = -.25f;
  float ratio = kRatioSentinel;
  const int original_workers = omp_get_max_threads();

  // The existing optimizer initializes converged=false. A zero iteration
  // budget tests that branch directly; production still requests100 iterations.
  const bool accepted = small_gicp::register_loop_gicp(
    source, target, identity, kVoxelMetres, kCovarianceNeighbors, kWorkerThreads,
    kMaxCorrespondenceMetres, 0, kHalfSourceRatio, output, ratio);

  CHECK_FALSE(accepted);
  CHECK(std::memcmp(&ratio, &kRatioSentinel, sizeof(float)) == 0);
  checkTransformBits(output, sentinel);
  CHECK(omp_get_max_threads() == original_workers);
  checkCloudUnchanged(*source, original_source);
  checkCloudUnchanged(*target, original_target);
}
