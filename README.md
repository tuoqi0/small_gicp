# small_gicp

C++17 point-cloud registration with ICP, point-to-plane ICP, GICP and VGICP.
The library provides native point storage, KD-tree and voxel maps, preprocessing,
registration factors, and serial/OpenMP implementations.

## Dependencies

- C++17 compiler and CMake 3.16 or newer.
- Installed the `robot_core` package. Its native math and algorithm targets are used directly.
- OpenMP for the production registration and self-fit targets.

## Build

Inside `slam_ws/src/small_gicp`, follow the
[workspace build instructions](../../README.md). The workspace builds the pinned
`src/robot_core` submodule and ROS packages together with
`colcon build --base-paths src`; no separate sibling checkout is needed.

For a standalone build, pass the installed `robot_core` prefix explicitly:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/path/to/robot_core-install -DCMAKE_INSTALL_PREFIX=/path/to/small_gicp-install
cmake --build build --parallel 2
cmake --install build
```

The package always builds the production library and `sc_self_fit` executable with OpenMP; no test or backend switch is required.

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

## Update Scan Context self-fit

`sc_self_fit` computes each stored keyframe's geometric fit at its recorded pose
and writes the scores into the existing Scan Context database. It uses native
Core PCD/SC codecs, prepared small_gicp clouds and Core geometric fit counts;
it does not optimize the poses. The installed executable is exported as the
CMake target `small_gicp::sc_self_fit`.

After the standalone install above, back up the database and run the CLI from
any working directory. A ROS or Python environment is not needed:

```sh
cp /path/to/map.pcd.sc /path/to/map.pcd.sc.backup
/path/to/small_gicp-install/bin/sc_self_fit \
  /path/to/map.pcd /path/to/map.pcd.sc /path/to/map.pcd.scans
```

The third path is optional and defaults to the second path plus `ans`, so
`map.pcd.sc` selects `map.pcd.scans`. The files must come from the same map and
contain the same ordered keyframes. Keel writes the companion clouds when
`scan_context.keep_scans` is enabled. The tool overwrites only the `.sc` file;
keep its backup because a failed write is not rolled back atomically. Map PCD
and companion scans remain inputs. Unusable keyframe clouds retain their
previous score; the printed processed count distinguishes them from scored
rows. Read/format errors and invalid leaf settings return exit1, usage returns
exit2, and successful saving returns exit0. Backend exceptions propagate.

`SF_LEAF` sets the voxel size in metres, with a default of `0.5`:

```sh
SF_LEAF=0.3 /path/to/small_gicp-install/bin/sc_self_fit \
  /path/to/map.pcd /path/to/map.pcd.sc
```

The existing numeric conversion accepts a numeric prefix; the resulting value
must be finite and positive. Geometry normals use a fixed `1.0 m` radius,
matching uses a fixed `1.0 m` distance and `30°` normal-angle gate, and source
normals with `|nz| > 0.7` are excluded. These settings do not read customized
node YAML. Scores describe the stored map/keyframe geometry and do not establish
relocalization accuracy for other online settings.

The native installation uses a relative runtime library path from its configured
binary directory to its configured library directory (`bin` to `lib` by default).
Keep those directories together when moving an installation. In a ROS workspace,
the relocalization package installs the same ELF at its libexec path:

```sh
source /opt/ros/jazzy/setup.bash
source /path/to/slam_ws/install/setup.bash
ros2 run small_gicp_relocalization sc_self_fit \
  /path/to/map.pcd /path/to/map.pcd.sc /path/to/map.pcd.scans
```

Source the workspace setup so the relocated ROS command finds the installed
small_gicp library. For manually separated prefixes, source ROS and the ROS
package setup, then prepend the actual native library directory, for example
`export LD_LIBRARY_PATH="/path/to/small_gicp-install/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"`.
Use the configured library directory if it differs from `lib`. The ROS copy
has the same contents as the native executable and is not a second build.

## 验证
生产配准通过 small_gicp_relocalization 与 Keel 的 ROS 链路验证；本包只构建并安装注册库和 sc_self_fit。
