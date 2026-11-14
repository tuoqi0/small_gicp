#pragma once

#include <alg/loop_registration_types.h>
#include <alg/point_cloud_types.h>

namespace small_gicp {

bool register_loop_gicp(
  const earth::mantle::alg::PointCloudXYZI::Ptr& src,
  const earth::mantle::alg::PointCloudXYZI::Ptr& tgt,
  const earth::mantle::alg::loop_closure::RegistrationTransform& guess,
  float voxel, int num_neighbors, int num_threads, float max_corr_dist,
  int max_iterations, float min_inlier_ratio,
  earth::mantle::alg::loop_closure::RegistrationTransform& out,
  float& inlier_ratio);

}  // namespace small_gicp
