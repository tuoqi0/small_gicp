# small_gicp

C++17 point-cloud registration with ICP, point-to-plane ICP, GICP and VGICP.
The library provides native point storage, KD-tree and voxel maps, preprocessing,
registration factors, and serial/OpenMP/TBB implementations.

## Dependencies

- C++17 compiler and CMake 3.16 or newer.
- Installed `slam_math` and `robot_core` packages.
- OpenMP for parallel helpers and tests.
- TBB when using the `*_tbb.hpp` interfaces.

## Build

Inside `slam_ws/src/small_gicp`, follow the
[workspace build instructions](../../README.md). The workspace builds the pinned
`src/robot_core` submodule and ROS packages together with
`colcon build --base-paths src`; no separate sibling checkout is needed.

For a standalone build, add the installed math packages to `CMAKE_PREFIX_PATH`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
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

## Test

```sh
ctest --test-dir build --output-on-failure
```

The suite covers registration, preprocessing, matrix systems, optimizer failures,
point/flat/voxel storage, ROS point serialization, projection and sorting.

To verify an installation from a separate build directory:

```sh
cmake -S test/install_consumer -B consumer-build \
  -DCMAKE_PREFIX_PATH=/path/to/install
cmake --build consumer-build --parallel 2
ctest --test-dir consumer-build --output-on-failure
```

The consumer reuses the public helper contract through the installed
`small_gicp::small_gicp` target.

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
