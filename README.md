# small_gicp

C++17 point-cloud registration with ICP, point-to-plane ICP, GICP and VGICP.
The library provides native point storage, KD-tree and voxel maps, preprocessing,
registration factors, and serial/OpenMP/TBB implementations.

## Dependencies

- C++17 compiler and CMake 3.16 or newer.
- Installed the `robot_core` package. Its native math and algorithm targets are used directly.
- OpenMP for parallel helpers and tests.
- TBB when using the `*_tbb.hpp` interfaces.

## Build

Inside `slam_ws/src/small_gicp`, follow the
[workspace build instructions](../../README.md). The workspace builds the pinned
`src/robot_core` submodule and ROS packages together with
`colcon build --base-paths src`; no separate sibling checkout is needed.

For a standalone build, pass the installed `robot_core` prefix explicitly:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/path/to/robot_core-install \
  -DCMAKE_INSTALL_PREFIX=/path/to/small_gicp-install
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
cmake --install build
```

The helper library is enabled by default. `BUILD_TESTING` controls native tests
and requires the test headers supplied by `robot_core`. Native TBB tests can be
built with `-DBUILD_WITH_TBB=ON` when TBB is installed.

## Use the native API

The helper accepts standard coordinate arrays. Transforms map source coordinates
into the target frame; matrix entries use `[row][column]`.

```cpp
#include <small_gicp/registration/registration_helper.hpp>
#include <math/fixed_matrix_ops.h>
#include <array>
#include <vector>

int main() {
  std::vector<std::array<double, 3>> target;
  for (int row = -5; row <= 5; ++row) {
    for (int col = -5; col <= 5; ++col) {
      target.push_back({{4, row * .25, col * .25}});
      target.push_back({{row * .25, 3, col * .25}});
      target.push_back({{row * .25, col * .25, -1}});
    }
  }
  auto source = target;
  for (auto& point : source) point[0] -= .08;

  small_gicp::RegistrationSetting setting;
  setting.num_threads = 2;
  setting.downsampling_resolution = .1;
  setting.max_correspondence_distance = 1.0;
  const auto initial = earth::mantle::math::identityMatrix<double, 4>();
  const auto result = small_gicp::align(target, source, initial, setting);
  // result.T_target_source: native 4x4 source-to-target transform.
  // result.H: native robot_core Hessian matrix (6x6).
  // result.num_inliers: number of accepted source points.
  return result.converged ? 0 : 1;
}
```

Save the example as `register_cloud.cpp` and link the installed target:

```cmake
cmake_minimum_required(VERSION 3.16)
project(register_cloud LANGUAGES CXX)
find_package(small_gicp CONFIG REQUIRED)
add_executable(register_cloud register_cloud.cpp)
target_link_libraries(register_cloud PRIVATE small_gicp::small_gicp)
target_compile_features(register_cloud PRIVATE cxx_std_17)
```

Set `CMAKE_PREFIX_PATH` to the installation prefix when configuring the consumer.
No source-tree include paths are needed.

For repeated registration, use `preprocess_points()` once and reuse the prepared
cloud and KD-tree. Lower-level `Registration<Factor, Reduction>` templates are
available for custom factors and correspondence searches. `PointCloud::point(i)`,
`normal(i)` and `cov(i)` return native storage references; traits return values.

For localization geometry verification, prepare KNN covariance and geometric
normals with a physical support radius:

```cpp
#include <small_gicp/registration/registration_preparation.hpp>
#include <array>
#include <vector>

std::vector<std::array<float, 3>> coordinates;
// Fill coordinates with XYZ in metres in one shared frame.
auto prepared = small_gicp::prepare_registration_cloud(
  coordinates, 0.25, 20, 4, 1.0);
```

The five arguments are XYZ coordinates, voxel size in metres, KNN neighbor
count, covariance/normal worker count, and geometry-normal radius in metres.
The OMP voxel sampler uses its own default of four workers. XYZ centroids and
all covariance entries retain float precision in native double buffers. The
radius normals face the coordinate origin; they are independent of the KNN
covariance used by GICP. The input coordinates remain unchanged, and empty input
returns a nonnull empty cloud. Unusable radius neighborhoods retain invalid
normals; the caller applies its existing geometric acceptance rules. An invalid
radius leaves the prepared covariance cloud's normal buffers untouched. Backend
exceptions propagate to the caller.

This API returns the prepared cloud; construct a KD-tree with the existing
native interfaces when needed. It does not choose a registration transform or
decide whether a localization result is acceptable. `preprocess_points()` uses
its own sampling and KNN-normal policy and remains useful for general alignment.

For loop-closure GICP, include
`<small_gicp/registration/loop_registration.hpp>` and call
`small_gicp::register_loop_gicp(source, target, guess, voxel, neighbors, threads,
max_correspondence_distance, max_iterations, minimum_inlier_ratio, output,
inlier_ratio)`. Source and target are `robot_core` native `PointCloudXYZI::Ptr`
clouds; the float 4×4 transforms use
`earth::mantle::alg::loop_closure::RegistrationTransform` from
`<alg/loop_registration_types.h>`. Link `small_gicp::small_gicp`; this entry point
does not require the loop-closure backend or ROS.

Voxel size and maximum correspondence distance are in metres. The sampler uses
the OMP path even with one worker, and the acceptance ratio uses the sampled
source count. Empty input or a nonconverged solve returns `false` without
changing either output. A converged solve writes `inlier_ratio` before checking
the requested minimum; if that check fails, the output transform remains
unchanged. On success the transform maps source to target. Inputs remain
unchanged and backend exceptions propagate.

## Test

```sh
ctest --test-dir build --output-on-failure
```

The suite covers registration, preprocessing, matrix systems, optimizer failures,
point/flat/voxel storage, ROS point serialization, projection and sorting.
`registration_preparation_native_test` checks the prepared cloud through actual
KD-tree lookup and GICP, float centroid/covariance storage, sparse and degenerate
neighborhoods, independent radius/KNN support, origin-facing normals, input
immutability and the existing invalid-radius output:

```sh
ctest --test-dir build --output-on-failure -R '^registration_preparation_native_test$'
```

`loop_registration_native_test` calls the loop GICP entry point with real native
clouds. It checks acceptance at exactly 50% support, rejection just above that
boundary, nonconvergence, output preservation, unchanged input payload and the
caller's OMP worker setting:

```sh
ctest --test-dir build --output-on-failure -R '^loop_registration_native_test$'
```

To verify an installation from a separate build directory:

```sh
cmake -S test/install_consumer -B consumer-build \
  -DCMAKE_PREFIX_PATH="/path/to/small_gicp-install;/path/to/robot_core-install"
cmake --build consumer-build --parallel 2
ctest --test-dir consumer-build --output-on-failure
```

The consumer runs `installed_native_helper`,
`installed_registration_preparation` and `installed_loop_registration` through
the installed `small_gicp::small_gicp` target. They reuse the native test contracts; the test
input is read from this checkout while all library headers and symbols come
from the selected installations.

## Documentation

See the [C++ documentation index](docs/index.md). The [paper](docs/paper.md) and
[historical benchmark results](BENCHMARK.md) describe the upstream publication.

## License
This package is released under the MIT license.

If you find this package useful for your project, please consider leaving a comment [here](https://github.com/koide3/small_gicp/issues/3). It would help the author receive recognition in his organization and keep working on this project.

Please also cite [the following article](https://joss.theoj.org/papers/10.21105/joss.06948) if you use this software in an academic work.

```
@article{small_gicp,
author = {Kenji Koide},
title = {{small\_gicp: Efficient and parallel algorithms for point cloud registration}},
journal = {Journal of Open Source Software},
month = aug,
number = {100},
pages = {6948},
volume = {9},
year = {2024},
doi = {10.21105/joss.06948}
}
```

## Contact

[Kenji Koide](https://staff.aist.go.jp/k.koide/), National Institute of Advanced Industrial Science and Technology (AIST)
