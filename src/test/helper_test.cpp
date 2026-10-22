#include <robot_core/test/doctest.h>
#include <string>
#include <cmath>
#include <small_gicp/registration/registration_helper.hpp>
#include <small_gicp/util/eigen.hpp>

#include <small_gicp/benchmark/read_points.hpp>
#include <small_gicp/registration/registration_helper.hpp>

using namespace small_gicp;

class HelperTest {
public:
  HelperTest() { setup(); }
  void setup() {
    // Load points
    target_raw = std::make_shared<PointCloud>(read_ply("data/target.ply"));
    source_raw = std::make_shared<PointCloud>(read_ply("data/source.ply"));
    for (const auto& point : target_raw->points) target_eigen.push_back(eigen_matrix(point));
    for (const auto& point : source_raw->points) source_eigen.push_back(eigen_matrix(point));

    std::ifstream ifs("data/T_target_source.txt");
    if (!ifs) {
      std::cerr << "error: failed to open T_target_source.txt" << std::endl;
    }
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        ifs >> T_target_source(i, j);
      }
    }
  }

  bool compare_transformation(const Eigen::Isometry3d& T1, const Eigen::Isometry3d& T2) {
    const Eigen::Isometry3d e = T1.inverse() * T2;
    const double error_rot = Eigen::AngleAxisd(e.linear()).angle();
    const double error_trans = e.translation().norm();

    const double rot_tol = 2.5 * M_PI / 180.0;
    const double trans_tol = 0.2;

    CHECK(std::abs(static_cast<double>(error_rot) - static_cast<double>(0.0)) <= (rot_tol));
    CHECK(std::abs(static_cast<double>(error_trans) - static_cast<double>(0.0)) <= (trans_tol));

    return error_rot < rot_tol && error_trans < trans_tol;
  }

protected:
  PointCloud::Ptr target_raw;         ///< Target points
  PointCloud::Ptr source_raw;         ///< Source points
  std::vector<Eigen::Vector4d> target_eigen;
  std::vector<Eigen::Vector4d> source_eigen;
  Eigen::Isometry3d T_target_source;  ///< Ground truth transformation
};

// Load check
TEST_CASE_FIXTURE(HelperTest, "HelperTest.LoadCheck") {
  CHECK(!(target_raw->empty()));
  CHECK(!(source_raw->empty()));
}

TEST_CASE_FIXTURE(HelperTest, "HelperTest.EmptyPreprocess") {
  auto empty_points = std::make_shared<PointCloud>();
  auto [target, target_tree] = preprocess_points(*empty_points, 0.1, 10, 1);
  CHECK((target));
  CHECK((target_tree));
  CHECK((target->size()) == (0));
}

TEST_CASE_FIXTURE(HelperTest, "HelperTest.Preprocess") {
  const std::vector<int> num_threads = {1, 2};

  for (int N : num_threads) {
    auto [target, target_tree] = preprocess_points(*target_raw, 0.1, 10, N);
    CHECK_MESSAGE((target), ("N="), (N));
    CHECK_MESSAGE((target_tree), ("N="), (N));
    CHECK_MESSAGE((target->size()) > (1000), ("N="), (N));
    CHECK_MESSAGE((target->size()) < (target_raw->size()), ("N="), (N));

    auto [target2, target_tree2] = preprocess_points(target_eigen, 0.1, 10, N);
    CHECK_MESSAGE((target2), ("N="), (N));
    CHECK_MESSAGE((target_tree2), ("N="), (N));
    CHECK_MESSAGE(std::abs(static_cast<double>(target->size()) - static_cast<double>(target2->size())) <= (100), ("N="), (N));
  }
}

TEST_CASE_FIXTURE(HelperTest, "HelperTest.EmptyGaussianVoxelMap") {
  auto empty_points = std::make_shared<PointCloud>();
  auto voxelmap = create_gaussian_voxelmap(*empty_points, 0.1);
  CHECK((voxelmap));
  CHECK((voxelmap->size()) == (0));
}

TEST_CASE_FIXTURE(HelperTest, "HelperTest.GaussianVoxelMap") {
  auto voxelmap = create_gaussian_voxelmap(*target_raw, 0.1);
  CHECK((voxelmap));
  CHECK((voxelmap->size()) > (1000));
  CHECK((voxelmap->size()) < (target_raw->size()));
}

TEST_CASE_FIXTURE(HelperTest, "HelperTest.AlignEigen") {
  for (const std::string parameter : {"ICP", "PLANE_ICP", "GICP", "VGICP"}) {
    SUBCASE(parameter.c_str()) {
      const std::string method = parameter;

      RegistrationSetting setting;
      if (method == "ICP") {
        setting.type = RegistrationSetting::ICP;
      } else if (method == "PLANE_ICP") {
        setting.type = RegistrationSetting::PLANE_ICP;
      } else if (method == "GICP") {
        setting.type = RegistrationSetting::GICP;
      } else if (method == "VGICP") {
        setting.type = RegistrationSetting::VGICP;
      } else {
        std::cerr << "error: unknown method " << method << std::endl;
      }

      // Forward test
      auto result = align(target_eigen, source_eigen, Eigen::Isometry3d::Identity(), setting);
      CHECK((compare_transformation(T_target_source, eigen_transform(result.T_target_source))));

      // Backward test
      auto result2 = align(source_eigen, target_eigen, Eigen::Isometry3d::Identity(), setting);
      CHECK((compare_transformation(T_target_source, eigen_transform(result2.T_target_source).inverse())));
    }
  }
}

TEST_CASE_FIXTURE(HelperTest, "HelperTest.AlignSmall") {
  for (const std::string parameter : {"ICP", "PLANE_ICP", "GICP", "VGICP"}) {
    SUBCASE(parameter.c_str()) {
      const std::string method = parameter;

      RegistrationSetting setting;
      if (method == "ICP") {
        setting.type = RegistrationSetting::ICP;
      } else if (method == "PLANE_ICP") {
        setting.type = RegistrationSetting::PLANE_ICP;
      } else if (method == "GICP") {
        setting.type = RegistrationSetting::GICP;
      } else if (method == "VGICP") {
        setting.type = RegistrationSetting::VGICP;
      } else {
        std::cerr << "error: unknown method " << method << std::endl;
      }

      RegistrationResult result;

      auto [target, target_tree] = preprocess_points(*target_raw, 0.1, 10, 1);
      auto [source, source_tree] = preprocess_points(*source_raw, 0.1, 10, 1);

      // Forward test
      if (method != "VGICP") {
        result = align(*target, *source, *target_tree, Eigen::Isometry3d::Identity(), setting);
      } else {
        auto target_voxelmap = create_gaussian_voxelmap(*target, 1.0);
        result = align(*target_voxelmap, *source, Eigen::Isometry3d::Identity(), setting);
      }
      CHECK((compare_transformation(T_target_source, eigen_transform(result.T_target_source))));

      // Backward test
      if (method != "VGICP") {
        result = align(*source, *target, *source_tree, Eigen::Isometry3d::Identity(), setting);
      } else {
        auto source_voxelmap = create_gaussian_voxelmap(*source, 1.0);
        result = align(*source_voxelmap, *target, Eigen::Isometry3d::Identity(), setting);
      }
      CHECK((compare_transformation(T_target_source, eigen_transform(result.T_target_source).inverse())));
    }
  }
}
