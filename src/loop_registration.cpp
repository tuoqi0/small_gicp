#include <small_gicp/registration/loop_registration.hpp>

#include <cstddef>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/ann/kdtree_omp.hpp>
#include <small_gicp/factors/gicp_factor.hpp>
#include <small_gicp/registration/registration.hpp>
#include <small_gicp/registration/reduction_omp.hpp>
#include <small_gicp/util/normal_estimation_omp.hpp>
#include <small_gicp/util/downsampling_omp.hpp>

namespace small_gicp {
using earth::mantle::alg::PointCloudXYZI;
using earth::mantle::alg::loop_closure::RegistrationTransform;

bool register_loop_gicp(
  const PointCloudXYZI::Ptr & src, const PointCloudXYZI::Ptr & tgt,
  const RegistrationTransform & guess, float voxel, int num_neighbors, int num_threads,
  float max_corr_dist, int max_iterations, float min_inlier_ratio,
  RegistrationTransform & out, float & inlier_ratio)
{
  // Use the inlier ratio rather than a standalone fitness score: repeated
  // geometry can make a large misalignment look numerically inexpensive.
  if (!src || !tgt || src->empty() || tgt->empty()) return false;

  auto toNative = [](const PointCloudXYZI & in) {
    small_gicp::PointCloud o;
    // Raw input needs only coordinates; covariance buffers are created for
    // downsampled points, not for every point in the full submap.
    o.points.resize(in.size());
    for (std::size_t i = 0; i < in.size(); ++i) {
      const auto & p = in[i];
      o.point(i) = small_gicp::PointCloud::NativePoint::from({p.x, p.y, p.z, 1.0});
    }
    return o;
  };

  // Keep the OMP sampling path even with one worker: preprocess_points() uses
  // a different serial sampler in that case and can change centroid ordering.
  auto s_cov = small_gicp::voxelgrid_sampling_omp(toNative(*src), voxel, num_threads);
  auto t_cov = small_gicp::voxelgrid_sampling_omp(toNative(*tgt), voxel, num_threads);
  if (s_cov->empty() || t_cov->empty()) return false;

  small_gicp::estimate_covariances_omp(*s_cov, num_neighbors, num_threads);
  small_gicp::estimate_covariances_omp(*t_cov, num_neighbors, num_threads);
  const small_gicp::KdTree<small_gicp::PointCloud> t_tree(
    t_cov, small_gicp::KdTreeBuilderOMP(num_threads));

  small_gicp::Registration<small_gicp::GICPFactor, small_gicp::ParallelReductionOMP> reg;
  reg.reduction.num_threads = num_threads;
  reg.rejector.max_dist_sq = max_corr_dist * max_corr_dist;
  reg.optimizer.max_iterations = max_iterations;

  small_gicp::Transform init;
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 4; ++col) init[row][col] = guess[row][col];
  const auto res = reg.align(*t_cov, *s_cov, t_tree, init);
  if (!res.converged) return false;

  const double ratio =
    s_cov->empty() ? 0.0 : static_cast<double>(res.num_inliers) / s_cov->size();
  inlier_ratio = static_cast<float>(ratio);
  if (ratio < min_inlier_ratio) return false;
  for (std::size_t row = 0; row < 4; ++row)
    for (std::size_t col = 0; col < 4; ++col)
      out[row][col] = static_cast<float>(res.T_target_source[row][col]);
  return true;
}

}  // namespace small_gicp
