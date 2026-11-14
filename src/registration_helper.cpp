// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#include <small_gicp/registration/registration_helper.hpp>

#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/ann/kdtree_omp.hpp>
#include <small_gicp/ann/gaussian_voxelmap.hpp>
#include <small_gicp/util/downsampling.hpp>
#include <small_gicp/util/downsampling_omp.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#include <small_gicp/util/normal_estimation_omp.hpp>

#include <small_gicp/factors/icp_factor.hpp>
#include <small_gicp/factors/plane_icp_factor.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/registration/reduction_omp.hpp>
#include <small_gicp/registration/registration.hpp>

namespace small_gicp {

// Preprocess points
std::pair<PointCloud::Ptr, std::shared_ptr<KdTree<PointCloud>>> preprocess_points(const PointCloud& points, double downsampling_resolution, int num_neighbors, int num_threads) {
  if (num_threads == 1) {
    auto downsampled = voxelgrid_sampling(points, downsampling_resolution);
    auto kdtree = std::make_shared<KdTree<PointCloud>>(downsampled);
    estimate_normals_covariances(*downsampled, *kdtree, num_neighbors);
    return {downsampled, kdtree};
  } else {
    auto downsampled = voxelgrid_sampling_omp(points, downsampling_resolution);
    auto kdtree = std::make_shared<KdTree<PointCloud>>(downsampled);
    estimate_normals_covariances_omp(*downsampled, *kdtree, num_neighbors, num_threads);
    return {downsampled, kdtree};
  }
}

// Create Gaussian voxel map
GaussianVoxelMap::Ptr create_gaussian_voxelmap(const PointCloud& points, double voxel_resolution) {
  auto voxelmap = std::make_shared<GaussianVoxelMap>(voxel_resolution);
  voxelmap->insert(points);
  return voxelmap;
}

// Align point clouds
RegistrationResult
align(const PointCloud& target, const PointCloud& source, const KdTree<PointCloud>& target_tree, const Transform& init_T, const RegistrationSetting& setting) {
  switch (setting.type) {
    default:
      std::cerr << "invalid registration type" << std::endl;
      abort();
    case RegistrationSetting::ICP: {
      Registration<ICPFactor, ParallelReductionOMP> registration;
      registration.reduction.num_threads = setting.num_threads;
      registration.rejector.max_dist_sq = setting.max_correspondence_distance * setting.max_correspondence_distance;
      registration.criteria.rotation_eps = setting.rotation_eps;
      registration.criteria.translation_eps = setting.translation_eps;
      registration.optimizer.max_iterations = setting.max_iterations;
      registration.optimizer.verbose = setting.verbose;
      return registration.align(target, source, target_tree, init_T);
    }
    case RegistrationSetting::PLANE_ICP: {
      Registration<PointToPlaneICPFactor, ParallelReductionOMP> registration;
      registration.reduction.num_threads = setting.num_threads;
      registration.rejector.max_dist_sq = setting.max_correspondence_distance * setting.max_correspondence_distance;
      registration.criteria.rotation_eps = setting.rotation_eps;
      registration.criteria.translation_eps = setting.translation_eps;
      registration.optimizer.max_iterations = setting.max_iterations;
      registration.optimizer.verbose = setting.verbose;
      return registration.align(target, source, target_tree, init_T);
    }
    case RegistrationSetting::GICP: {
      Registration<GICPFactor, ParallelReductionOMP> registration;
      registration.reduction.num_threads = setting.num_threads;
      registration.rejector.max_dist_sq = setting.max_correspondence_distance * setting.max_correspondence_distance;
      registration.criteria.rotation_eps = setting.rotation_eps;
      registration.criteria.translation_eps = setting.translation_eps;
      registration.optimizer.max_iterations = setting.max_iterations;
      registration.optimizer.verbose = setting.verbose;
      return registration.align(target, source, target_tree, init_T);
    }
    case RegistrationSetting::VGICP: {
      std::cerr << "error: use align(const GaussianVoxelMap&, const GaussianVoxelMap&, const Transform&, const RegistrationSetting&) for VGICP" << std::endl;
      return RegistrationResult{};
    }
  }
}

// Align point clouds with VGICP
RegistrationResult align(const GaussianVoxelMap& target, const PointCloud& source, const Transform& init_T, const RegistrationSetting& setting) {
  if (setting.type != RegistrationSetting::VGICP) {
    std::cerr << "invalid registration type for GaussianVoxelMap" << std::endl;
  }

  Registration<GICPFactor, ParallelReductionOMP> registration;
  registration.reduction.num_threads = setting.num_threads;
  registration.criteria.rotation_eps = setting.rotation_eps;
  registration.criteria.translation_eps = setting.translation_eps;
  registration.optimizer.max_iterations = setting.max_iterations;
  registration.optimizer.verbose = setting.verbose;
  return registration.align(target, source, target, init_T);
}

}  // namespace small_gicp
