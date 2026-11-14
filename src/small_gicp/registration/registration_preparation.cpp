#include <small_gicp/registration/registration_preparation.hpp>

#include <alg/fpfh.h>
#include <alg/point_cloud_types.h>
#include <cstddef>
#include <small_gicp/util/downsampling_omp.hpp>
#include <small_gicp/util/normal_estimation_omp.hpp>

namespace small_gicp {

PointCloud::Ptr prepare_registration_cloud(
  const std::vector<std::array<float, 3>>& coordinates, double leaf,
  int neighbors, int threads, double normal_radius)
{
  const PointCloud input(coordinates);
  auto cloud = voxelgrid_sampling_omp(input, leaf);
  for (auto& point : cloud->points)
    for (std::size_t axis = 0; axis < 3; ++axis)
      point[axis][0] = static_cast<float>(point[axis][0]);
  estimate_covariances_omp(*cloud, neighbors, threads);
  for (auto& covariance : cloud->covs)
    for (std::size_t i = 0; i < 16; ++i)
      covariance.raw_pointer()[i] = static_cast<float>(covariance.raw_pointer()[i]);
  earth::mantle::alg::PointCloudXYZI normal_input;
  normal_input.reserve(cloud->size());
  for (const auto& point : cloud->points)
    normal_input.push_back(
      {static_cast<float>(point[0][0]), static_cast<float>(point[1][0]),
       static_cast<float>(point[2][0]), 0});
  const auto geometry =
    earth::mantle::alg::estimateRadiusNormals(normal_input, normal_radius, {{0, 0, 0}}, threads);
  for (std::size_t i = 0; i < geometry.size(); ++i) {
    const auto& normal = geometry[i];
    cloud->normal(i) = PointCloud::NativePoint::from(
      {normal.normal_x, normal.normal_y, normal.normal_z, 0});
  }
  return cloud;
}

}  // namespace small_gicp
