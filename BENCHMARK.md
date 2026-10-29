# Benchmark

The results and plots below are retained upstream benchmarks. They were not
re-measured for this fork's current native backend and do not establish its
performance or accuracy. For the native workspace build and validation entry,
see the [workspace instructions](../../README.md) and this fork's
[dependency requirements](README.md#dependencies).

## Dataset

The following benchmark requires the KITTI odometry evaluation dataset. You can download the full dataset (80GB) from [the official dataset page](https://www.cvlibs.net/datasets/kitti/eval_odometry.php) or a part of the dataset (500 frames in 00 sequence, 622MB) from [google drive (KITTI00.tar.gz)](https://drive.google.com/file/d/1h9tARKvX6BwLfc_vfMfdxuP3bSbmgjmd/view?usp=sharing).

Note that because the original KITTI dataset is distributed under the CC BY-NC-SA 3.0 license, the derived dataset (KITTI00.tar.gz) must not be used for commercial purposes.

## Current native delivery

The upstream benchmark executables, comparison adapters and launcher scripts are
not included in this native-only delivery. Their former build commands do not
apply. Build and validate the native library using the [README](README.md).
The following results and assets remain historical publication material.

## Results

All benchmarks were conducted on the KITTI 00 sequence.

### Downsampling

- Single-threaded `small_gicp::voxelgrid_sampling` is about **1.3x faster** than `pcl::VoxelGrid`.
- Multi-threaded `small_gicp::voxelgrid_sampling_tbb` (6 threads) is about **3.2x faster** than `pcl::VoxelGrid`.
- `small_gicp::voxelgrid_sampling` gives accurate downsampling results (almost identical to those of `pcl::VoxelGrid`) while `pcl::ApproximateVoxelGrid` yields spurious points (up to 2x points).
- `small_gicp::voxelgrid_sampling` can process a larger point cloud with a fine voxel resolution compared to `pcl::VoxelGrid` (for a point cloud of 1000m width, minimum voxel resolution can be 0.5 mm).

![downsampling_comp](docs/assets/downsampling_comp.png)

- While TBB shows slightly better scalability, both the parallelism backends do not obtain a speed gain for the cases with threads more than 16.

![downsampling_threads](docs/assets/downsampling_threads.png)

### KdTree construction

- Multi-threaded implementation (TBB and OMP) can be up to **4x faster** than the single-threaded one (All the implementations are based on nanoflann).
- ~~The processing speed gets faster as the number of threads increases, but the speed gain is not monotonic sometimes (because of the scheduling algorithm or some CPU(AMD 5995WX)-specific issues?)~~.
- The new KdTree implementation shows a good scalability thanks to its well balanced task assignment.
- This benchmark only compares the construction time (query time is not included). 

![kdtree_time](docs/assets/kdtree_time.png)

### Odometry estimation

- Single-thread `small_gicp::GICP` is about **2.4x and 1.9x faster** than `pcl::GICP` and `fast_gicp::GICP`, respectively.
- `small_gicp::(GICP|VGICP)` shows a better multi-thread scalability compared to `fast_gicp::(GICP|VGICP)`.
- `small_gicp::GICP` parallelized with TBB flow graph implementation shows an excellent scalability to many-threads situations (**~128 threads**) but with latency degradation.

![odometry_time](docs/assets/odometry_time.png)

**SIMD intrinsics (-march=native)** (We recommend keeping this feature disabled unless you are 100% sure what it is)

- `BUILD_WITH_MARCH_NATIVE=ON` enables platform-specific intrinsics and squeezing the performance (**1.1x speedup for free**).
- However, you must ensure that all involved libraries are built with `-march=native`, otherwise the program will crash.
- Generally, it is difficult to properly set `-march=native` for all libraries, and we recommend keeping `BUILD_WITH_MARCH_NATIVE=OFF`.

Results:
- `BUILD_WITH_MARCH_NATIVE=OFF` : `Eigen::SimdInstructionSetsInUse()=SSE, SSE2`
- `BUILD_WITH_MARCH_NATIVE=ON` : `Eigen::SimdInstructionSetsInUse()=AVX, SSE, SSE2, SSE3, SSSE3, SSE4.1, SSE4.2`

![odometry_native](docs/assets/odometry_native.png)

**Accuracy**

- `small_gicp::GICP` outputs mostly identical results to those of `fast_gicp::GICP`.
- The results of `small_gicp::VGICP` slightly differ from `fast_gicp::VGICP`. Although the difference is marginal, it needs to be investigated.

```
pcl_gicp             : APE=6.451 +- 3.421  RPE(100)=2.424 +- 1.707  RPE(400)=8.416 +- 4.284  RPE(800)=12.652 +- 6.799
fast_gicp            : APE=6.118 +- 3.078  RPE(100)=1.212 +- 0.717  RPE(400)=6.058 +- 3.128  RPE(800)=10.356 +- 6.335
fast_vgicp           : APE=6.791 +- 3.215  RPE(100)=1.253 +- 0.734  RPE(400)=6.315 +- 3.011  RPE(800)=10.367 +- 6.147
small_gicp           : APE=6.096 +- 3.056  RPE(100)=1.211 +- 0.717  RPE(400)=6.057 +- 3.123  RPE(800)=10.364 +- 6.336
small_gicp (tbb)     : APE=6.096 +- 3.056  RPE(100)=1.211 +- 0.717  RPE(400)=6.057 +- 3.123  RPE(800)=10.364 +- 6.336
small_gicp (omp)     : APE=6.096 +- 3.056  RPE(100)=1.211 +- 0.717  RPE(400)=6.057 +- 3.123  RPE(800)=10.364 +- 6.336
small_vgicp          : APE=5.956 +- 2.725  RPE(100)=1.315 +- 0.762  RPE(400)=6.849 +- 3.401  RPE(800)=10.396 +- 6.972
```

### Comparison with Open3D

[Code](https://github.com/koide3/small_gicp/blob/pybench/src/benchmark/odometry_benchmark.py)

Processing speed comparison between small_gicp and Open3D ([youtube](https://youtu.be/LNESzGXPr4c?feature=shared)).
[![small_comp](https://github.com/koide3/small_gicp/assets/31344317/7959edd6-f0e4-4318-b4c1-a3f8755c407f)](https://youtu.be/LNESzGXPr4c?feature=shared)
