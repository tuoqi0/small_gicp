#include <small_gicp/util/eigen.hpp>
#include <robot_core/test/doctest.h>
#include <string>
#include <cmath>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/features/normal_3d_omp.h>

#include <small_gicp/pcl/pcl_point.hpp>
#include <small_gicp/pcl/pcl_point_traits.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/ann/kdtree.hpp>
#include <small_gicp/util/downsampling.hpp>
#include <small_gicp/util/normal_estimation.hpp>
#include <small_gicp/util/normal_estimation_omp.hpp>
#include <small_gicp/util/normal_estimation_tbb.hpp>
#include <small_gicp/benchmark/read_points.hpp>

using namespace small_gicp;

class NormalEstimationTest {
public:
  NormalEstimationTest() { setup(); }
  void setup() {
    points = std::make_shared<PointCloud>(read_ply("data/target.ply"));
    points = voxelgrid_sampling(*points, 0.25);
    estimate_normals_covariances(*points, num_neighbors);

    tree = std::make_shared<KdTree<PointCloud>>(points);

    points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointNormal>>();
    points_pcl->resize(points->size());
    for (size_t i = 0; i < points->size(); i++) {
      points_pcl->at(i).getVector4fMap() = small_gicp::eigen_view(points->point(i)).cast<float>();
    }

    pcl::NormalEstimationOMP<pcl::PointNormal, pcl::PointNormal> normal_estimation;
    normal_estimation.setKSearch(num_neighbors);
    normal_estimation.setInputCloud(points_pcl);
    normal_estimation.compute(*points_pcl);
  }

  template <typename PointCloud>
  bool check_normals(const PointCloud& points) {
    if (traits::size(points) != points_pcl->size()) {
      return false;
    }

    for (size_t i = 0; i < traits::size(points); i++) {
      const Eigen::Vector4d normal = small_gicp::eigen_matrix(traits::normal(points, i));
      CHECK_MESSAGE(std::abs(static_cast<double>(normal.w()) - static_cast<double>(0.0)) <= (1e-6), ("Last element must be zero"));
      CHECK_MESSAGE(std::abs(static_cast<double>(normal.norm()) - static_cast<double>(1.0)) <= (1e-6), ("Normal must be unit vector"));

      const double error = (normal - points_pcl->at(i).getNormalVector4fMap().cast<double>()).norm();
      CHECK_MESSAGE(std::abs(static_cast<double>(error) - static_cast<double>(0.0)) <= (1e-3), ("Normal estimation error is too large"));
    }
    return true;
  }

  template <typename PointCloud>
  bool check_covs(const PointCloud& points) {
    if (traits::size(points) != this->points->size()) {
      return false;
    }

    for (size_t i = 0; i < traits::size(points); i++) {
      const Eigen::Matrix4d cov = small_gicp::eigen_matrix(traits::cov(points, i));
      CHECK_MESSAGE(std::abs(static_cast<double>(cov.bottomRows(1).norm()) - static_cast<double>(0.0)) <= (1e-6), ("Bottom row must be zero"));
      CHECK_MESSAGE(std::abs(static_cast<double>(cov.rightCols(1).norm()) - static_cast<double>(0.0)) <= (1e-6), ("Right column must be zero"));

      for (int j = 0; j < 3; j++) {
        for (int k = j + 1; k < 3; k++) {
          CHECK_MESSAGE(std::abs(static_cast<double>(cov(j, k)) - static_cast<double>(cov(k, j))) <= (1e-6), ("Covariance matrix must be symmetric"));
        }
      }

      const double error = (cov - small_gicp::eigen_view(this->points->cov(i))).norm();
      CHECK_MESSAGE(std::abs(static_cast<double>(error) - static_cast<double>(0.0)) <= (1e-3), ("Covariance estimation error is too large"));
    }
    return true;
  }

protected:
  const int num_neighbors = 20;

  PointCloud::Ptr points;
  KdTree<PointCloud>::Ptr tree;
  pcl::PointCloud<pcl::PointNormal>::Ptr points_pcl;
};

// Load check
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.LoadCheck") {
  CHECK(!(points->empty()));
  CHECK(!(points_pcl->empty()));
}

// Empty test
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.EmptyTest") {
  // Empty point cloud
  auto empty_points = std::make_shared<PointCloud>();
  auto empty_tree = std::make_shared<KdTree<PointCloud>>(empty_points);

  estimate_normals(*empty_points, *empty_tree, num_neighbors);
  CHECK((empty_points->empty()));
  estimate_covariances(*empty_points, *empty_tree, num_neighbors);
  CHECK((empty_points->empty()));
  estimate_normals_covariances(*empty_points, *empty_tree, num_neighbors);
  CHECK((empty_points->empty()));

  // Empty point cloud (PCL)
  auto empty_points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointNormalCovariance>>();
  auto empty_tree_pcl = std::make_shared<KdTree<pcl::PointCloud<pcl::PointNormalCovariance>>>(empty_points_pcl);

  estimate_normals(*empty_points_pcl, *empty_tree_pcl, num_neighbors);
  CHECK((empty_points_pcl->empty()));
  estimate_covariances(*empty_points_pcl, *empty_tree_pcl, num_neighbors);
  CHECK((empty_points_pcl->empty()));
  estimate_normals_covariances(*empty_points_pcl, *empty_tree_pcl, num_neighbors);
  CHECK((empty_points_pcl->empty()));
}

// Empty test (TBB)
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.EmptyTestTBB") {
  // Empty point cloud
  auto empty_points = std::make_shared<PointCloud>();
  auto empty_tree = std::make_shared<KdTree<PointCloud>>(empty_points);

  estimate_normals_tbb(*empty_points, *empty_tree, num_neighbors);
  CHECK((empty_points->empty()));
  estimate_covariances_tbb(*empty_points, *empty_tree, num_neighbors);
  CHECK((empty_points->empty()));
  estimate_normals_covariances_tbb(*empty_points, *empty_tree, num_neighbors);
  CHECK((empty_points->empty()));

  // Empty point cloud (PCL)
  auto empty_points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointNormalCovariance>>();
  auto empty_tree_pcl = std::make_shared<KdTree<pcl::PointCloud<pcl::PointNormalCovariance>>>(empty_points_pcl);

  estimate_normals_tbb(*empty_points_pcl, *empty_tree_pcl, num_neighbors);
  CHECK((empty_points_pcl->empty()));
  estimate_covariances_tbb(*empty_points_pcl, *empty_tree_pcl, num_neighbors);
  CHECK((empty_points_pcl->empty()));
  estimate_normals_covariances_tbb(*empty_points_pcl, *empty_tree_pcl, num_neighbors);
  CHECK((empty_points_pcl->empty()));
}

// Empty test (OMP)
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.EmptyTestOMP") {
  // Empty point cloud
  auto empty_points = std::make_shared<PointCloud>();
  auto empty_tree = std::make_shared<KdTree<PointCloud>>(empty_points);

  estimate_normals_omp(*empty_points, *empty_tree, num_neighbors, 2);
  CHECK((empty_points->empty()));
  estimate_covariances_omp(*empty_points, *empty_tree, num_neighbors, 2);
  CHECK((empty_points->empty()));
  estimate_normals_covariances_omp(*empty_points, *empty_tree, num_neighbors, 2);
  CHECK((empty_points->empty()));

  // Empty point cloud (PCL)
  auto empty_points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointNormalCovariance>>();
  auto empty_tree_pcl = std::make_shared<KdTree<pcl::PointCloud<pcl::PointNormalCovariance>>>(empty_points_pcl);

  estimate_normals_omp(*empty_points_pcl, *empty_tree_pcl, num_neighbors, 2);
  CHECK((empty_points_pcl->empty()));
  estimate_covariances_omp(*empty_points_pcl, *empty_tree_pcl, num_neighbors, 2);
  CHECK((empty_points_pcl->empty()));
  estimate_normals_covariances_omp(*empty_points_pcl, *empty_tree_pcl, num_neighbors, 2);
  CHECK((empty_points_pcl->empty()));
}

// Normal/covariance estimation test
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.NormalEstimationTest") {
  auto estimated = std::make_shared<PointCloud>();
  *estimated = *points;

  estimate_normals(*estimated, *tree, num_neighbors);
  CHECK((check_normals(*estimated)));

  estimate_covariances(*estimated, *tree, num_neighbors);
  CHECK((check_covs(*estimated)));

  *estimated = *points;
  estimate_normals_covariances(*estimated, *tree, num_neighbors);
  CHECK((check_normals(*estimated)));
  CHECK((check_covs(*estimated)));
}

// Normal/covariance estimation test (TBB)
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.NormalEstimationTestTBB") {
  auto estimated = std::make_shared<PointCloud>();
  *estimated = *points;

  estimate_normals_tbb(*estimated, *tree, num_neighbors);
  CHECK((check_normals(*estimated)));

  estimate_covariances_tbb(*estimated, *tree, num_neighbors);
  CHECK((check_covs(*estimated)));

  *estimated = *points;
  estimate_normals_covariances_tbb(*estimated, *tree, num_neighbors);
  CHECK((check_normals(*estimated)));
  CHECK((check_covs(*estimated)));
}

// Normal/covariance estimation test (OMP)
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.NormalEstimationTestOMP") {
  auto estimated = std::make_shared<PointCloud>();
  *estimated = *points;

  estimate_normals_omp(*estimated, *tree, num_neighbors, 2);
  CHECK((check_normals(*estimated)));

  estimate_covariances_omp(*estimated, *tree, num_neighbors, 2);
  CHECK((check_covs(*estimated)));

  *estimated = *points;
  estimate_normals_covariances_omp(*estimated, *tree, num_neighbors, 2);
  CHECK((check_normals(*estimated)));
  CHECK((check_covs(*estimated)));
}

// Normal/covariance estimation test (PCL)
TEST_CASE_FIXTURE(NormalEstimationTest, "NormalEstimationTest.NormalEstimationTestPCL") {
  const auto copy_points = [&](auto& points) {
    points.resize(points_pcl->size());
    for (size_t i = 0; i < points_pcl->size(); i++) {
      points.at(i).getVector4fMap() = points_pcl->at(i).getVector4fMap();
    }
  };

  // Normal estimation
  auto point_normals = pcl::make_shared<pcl::PointCloud<pcl::PointNormal>>();
  copy_points(*point_normals);

  estimate_normals(*point_normals, num_neighbors);
  CHECK((check_normals(*point_normals)));
  estimate_normals_tbb(*point_normals, num_neighbors);
  CHECK((check_normals(*point_normals)));
  estimate_normals_omp(*point_normals, num_neighbors, 2);
  CHECK((check_normals(*point_normals)));

  // Covariance estimation
  auto point_covs = pcl::make_shared<pcl::PointCloud<pcl::PointCovariance>>();
  copy_points(*point_covs);

  estimate_covariances(*point_covs, num_neighbors);
  CHECK((check_covs(*point_covs)));
  estimate_covariances_tbb(*point_covs, num_neighbors);
  CHECK((check_covs(*point_covs)));
  estimate_covariances_omp(*point_covs, num_neighbors, 2);
  CHECK((check_covs(*point_covs)));

  // Normal/covariance estimation
  auto point_normals_covs = pcl::make_shared<pcl::PointCloud<pcl::PointNormalCovariance>>();
  copy_points(*point_normals_covs);

  estimate_normals_covariances(*point_normals_covs, num_neighbors);
  CHECK((check_normals(*point_normals_covs)));
  CHECK((check_covs(*point_normals_covs)));
  estimate_normals_covariances_tbb(*point_normals_covs, num_neighbors);
  CHECK((check_normals(*point_normals_covs)));
  CHECK((check_covs(*point_normals_covs)));
  estimate_normals_covariances_omp(*point_normals_covs, num_neighbors, 2);
  CHECK((check_normals(*point_normals_covs)));
  CHECK((check_covs(*point_normals_covs)));
}
