#include <unordered_set>

#include <robot_core/test/doctest.h>
#include <small_gicp/util/eigen.hpp>
#include <string>
#include <cmath>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/voxel_grid.h>

#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/pcl/pcl_point.hpp>
#include <small_gicp/pcl/pcl_point_traits.hpp>
#include <small_gicp/util/downsampling.hpp>
#include <small_gicp/util/downsampling_omp.hpp>
#include <small_gicp/util/downsampling_tbb.hpp>
#include <small_gicp/benchmark/read_points.hpp>

using namespace small_gicp;

class DownsamplingTest {
public:
  std::string selected_method_;
  DownsamplingTest() { setup(); }
  void setup() {
    // Load points
    auto points_4f = read_ply("data/target.ply");
    points = std::make_shared<PointCloud>(points_4f);
    points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
    points_pcl->resize(points_4f.size());
    for (size_t i = 0; i < points_4f.size(); i++) {
      points_pcl->at(i).getVector4fMap() = points_4f[i];
    }

    // Downsample points using pcl::VoxelGrid for reference
    resolutions = {0.1, 0.5, 1.0};
    for (double resolution : resolutions) {
      pcl::VoxelGrid<pcl::PointXYZ> voxelgrid;
      voxelgrid.setLeafSize(resolution, resolution, resolution);
      voxelgrid.setInputCloud(points_pcl);

      auto downsampled = pcl::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
      voxelgrid.filter(*downsampled);
      downsampled_pcl.push_back(downsampled);
    }
  }

  // Apply downsampling
  template <typename PointCloud>
  std::shared_ptr<PointCloud> downsample(const PointCloud& points, double resolution) {
    const std::string method = selected_method_;
    if (method == "SMALL") {
      return voxelgrid_sampling(points, resolution);
    } else if (method == "TBB") {
      return voxelgrid_sampling_tbb(points, resolution);
    } else if (method == "OMP") {
      return voxelgrid_sampling_omp(points, resolution);
    } else {
      throw std::runtime_error("Invalid method: " + method);
    }
  }

protected:
  std::vector<double> resolutions;                                   ///< Downsampling resolutions
  PointCloud::Ptr points;                                            ///< Input points
  pcl::PointCloud<pcl::PointXYZ>::Ptr points_pcl;                    ///< Input points (pcl)
  std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> downsampled_pcl;  ///< Reference downsampling results (pcl)
};

// Load check
TEST_CASE_FIXTURE(DownsamplingTest, "DownsamplingTest.LoadCheck") {
  CHECK(!(points->empty()));
  CHECK(!(points_pcl->empty()));
  for (auto downsampled : downsampled_pcl) {
    CHECK(!(downsampled->empty()));
  }
}

// Check if downsampling works correctly for empty points
TEST_CASE_FIXTURE(DownsamplingTest, "DownsamplingTest.EmptyTest") {
  for (const std::string parameter : {"SMALL", "TBB", "OMP"}) {
    SUBCASE(parameter.c_str()) {
      selected_method_ = parameter;
      auto empty_points = std::make_shared<PointCloud>();
      auto empty_downsampled = downsample(*empty_points, 0.1);
      CHECK((empty_downsampled));
      CHECK_MESSAGE((empty_downsampled->size()) == (0), ("Empty test small: " + parameter));

      auto empty_points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
      auto empty_downsampled_pcl = downsample(*empty_points_pcl, 0.1);
      CHECK((empty_downsampled_pcl));
      CHECK_MESSAGE((empty_downsampled_pcl->size()) == (0), ("Empty test pcl: " + parameter));
    }
  }
}

// Check if downsampling results are mostly identical to those of pcl::VoxelGrid
TEST_CASE_FIXTURE(DownsamplingTest, "DownsamplingTest.DownsampleTest") {
  for (const std::string parameter : {"SMALL", "TBB", "OMP"}) {
    SUBCASE(parameter.c_str()) {
      selected_method_ = parameter;
      for (size_t i = 0; i < resolutions.size(); i++) {
        auto result = downsample(*points, resolutions[i]);
        CHECK_MESSAGE((std::abs(1.0 - static_cast<double>(result->size()) / downsampled_pcl[i]->size())) < (0.9), ("Downsampled size check (small): " + parameter));
        auto result_pcl = downsample(*points_pcl, resolutions[i]);
        CHECK_MESSAGE((std::abs(1.0 - static_cast<double>(result_pcl->size()) / downsampled_pcl[i]->size())) < (0.9), ("Downsampled size check (pcl): " + parameter));
        CHECK_MESSAGE((result->size()) == (result_pcl->size()), ("Size check (small vs pcl): " + parameter));
      }
    }
  }
}

// Check if random sampling works correctly for empty points
TEST_CASE_FIXTURE(DownsamplingTest, "DownsamplingTest.EmptyRandamSamplingTest") {
  for (const std::string parameter : {"SMALL", "TBB", "OMP"}) {
    SUBCASE(parameter.c_str()) {
      selected_method_ = parameter;
      std::mt19937 mt;

      auto empty_points = std::make_shared<PointCloud>();
      auto empty_downsampled = random_sampling(*empty_points, 1000, mt);
      CHECK((empty_downsampled));
      CHECK_MESSAGE((empty_downsampled->size()) == (0), ("Empty test small: " + parameter));

      auto empty_points_pcl = pcl::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
      auto empty_downsampled_pcl = random_sampling(*empty_points_pcl, 1000, mt);
      CHECK((empty_downsampled_pcl));
      CHECK_MESSAGE((empty_downsampled_pcl->size()) == (0), ("Empty test pcl: " + parameter));
    }
  }
}

// Test random sampling
TEST_CASE_FIXTURE(DownsamplingTest, "DownsamplingTest.RandamSamplingTest") {
  for (const std::string parameter : {"SMALL", "TBB", "OMP"}) {
    SUBCASE(parameter.c_str()) {
      selected_method_ = parameter;
      std::mt19937 mt;

      auto downsampled = voxelgrid_sampling(*points, 0.1);

      const std::vector<size_t> num_points = {0, 100, 1000};
      for (size_t N : num_points) {
        auto result = random_sampling(*downsampled, N, mt);
        CHECK((result));
        CHECK_MESSAGE((result->size()) == (N), ("Size check (small): " + parameter));

        std::unordered_set<size_t> indices;
        for (size_t i = 0; i < N; i++) {
          const auto found = std::find_if(downsampled->points.begin(), downsampled->points.end(), [&](const auto& p) {
            return (small_gicp::eigen_view(p) - small_gicp::eigen_view(result->points[i])).norm() < 1.0e-6;
          });
          CHECK_MESSAGE((found) != (downsampled->points.end()), ("Existence check (small): " + parameter));

          const size_t index = std::distance(downsampled->points.begin(), found);
          CHECK_MESSAGE((indices.count(index)) == (0), ("Uniqueness check (small): " + parameter));
          indices.insert(index);
        }
      }
    }
  }
}
