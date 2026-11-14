// Computes per-keyframe self-fit for a scan-context database and writes it back.
//
// Self-fit is how well a keyframe's own cloud matches the map at that keyframe's
// own recorded pose: it measures what the place offers, not whether a match is
// correct.
//
// Usage: sc_self_fit <map.pcd> <map.pcd.sc> [map.pcd.scans]
// The scans file defaults to <sc path> + "ans" and is written by keel_lio when
// scan_context.keep_scans is on.

#include <alg/point_cloud_io.h>
#include <alg/registration_verification.h>
#include <alg/scan_context.h>
#include <math/fixed_matrix_ops.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <small_gicp/ann/kdtree_omp.hpp>
#include <small_gicp/points/point_cloud.hpp>
#include <small_gicp/registration/registration_preparation.hpp>
#include <string>
#include <vector>

namespace scan_context = earth::mantle::alg::scan_context;

namespace
{

using Cloud = small_gicp::PointCloud;
// Uses the node's default feature_normal_radius in metres; this offline tool
// does not read customized node parameters.
constexpr double kNormalRadiusMeters = 1.0;

}  // namespace

int main(int argc, char ** argv)
{
  if (argc < 3) {
    std::fprintf(stderr, "usage: sc_self_fit <map.pcd> <map.pcd.sc> [map.pcd.scans]\n");
    return 2;
  }
  const std::string map_path = argv[1];
  const std::string sc_path = argv[2];
  const std::string scans_path = argc > 3 ? argv[3] : sc_path + "ans";

  scan_context::Database db;
  if (!db.load(sc_path)) {
    std::fprintf(stderr, "cannot read scan-context database: %s\n", sc_path.c_str());
    return 1;
  }
  std::vector<std::vector<float>> scans;
  if (!scan_context::loadScans(scans_path, scans)) {
    std::fprintf(stderr, "cannot read keyframe clouds: %s\n", scans_path.c_str());
    return 1;
  }
  if (scans.size() != db.entries.size()) {
    std::fprintf(
      stderr, "entry count mismatch: %zu clouds vs %zu database entries\n", scans.size(),
      db.entries.size());
    return 1;
  }

  earth::mantle::alg::PointCloud<earth::mantle::alg::PointXYZI> raw;
  if (earth::mantle::alg::loadPCDFile(map_path, raw) < 0 || raw.empty()) {
    std::fprintf(stderr, "cannot read map: %s\n", map_path.c_str());
    return 1;
  }
  const double leaf = std::getenv("SF_LEAF") ? std::atof(std::getenv("SF_LEAF")) : 0.5;
  if (!std::isfinite(leaf) || leaf <= 0) {
    std::fprintf(stderr, "SF_LEAF must be finite and positive\n");
    return 1;
  }
  std::vector<std::array<float, 3>> coordinates;
  coordinates.reserve(raw.size());
  for (const auto & point : raw) coordinates.push_back({{point.x, point.y, point.z}});
  auto target =
    small_gicp::prepare_registration_cloud(coordinates, leaf, 10, 8, kNormalRadiusMeters);
  auto tree = std::make_shared<small_gicp::KdTree<Cloud>>(target, small_gicp::KdTreeBuilderOMP(8));
  std::printf("map %zu points, %zu keyframes\n", target->size(), db.entries.size());

  const double d2 = 1.0;  // This offline tool uses a 1 m distance gate.
  const double cos_min = std::cos(30.0 * M_PI / 180.0);
  std::size_t done = 0;
  for (std::size_t i = 0; i < db.entries.size(); ++i) {
    const auto & sv = scans[i];
    if (sv.size() < 3) continue;
    std::vector<std::array<float, 3>> source;
    source.reserve(sv.size() / 3);
    // Clouds are stored recentred on the entry; put them back in the map frame.
    for (std::size_t k = 0; k + 2 < sv.size(); k += 3) {
      source.push_back(
        {{sv[k] + db.entries[i].x, sv[k + 1] + db.entries[i].y, sv[k + 2] + db.entries[i].z}});
    }
    auto src = small_gicp::prepare_registration_cloud(source, leaf, 10, 8, kNormalRadiusMeters);
    if (src->size() < 20) continue;

    const earth::mantle::alg::RegistrationFitCriteria criteria{d2, 0.7, cos_min, true};
    const auto counts = earth::mantle::alg::registrationFitCounts(
      *src, target.get(), *tree, earth::mantle::math::identityMatrix<double, 4>(), criteria);
    db.entries[i].self_fit =
      counts.considered
        ? static_cast<float>(
            static_cast<double>(counts.matched) / static_cast<double>(counts.considered))
        : 0.0f;
    ++done;
  }

  if (!db.save(sc_path)) {
    std::fprintf(stderr, "write back failed: %s\n", sc_path.c_str());
    return 1;
  }
  double sum = 0.0, lo = 1.0, hi = 0.0;
  for (const auto & e : db.entries) {
    sum += e.self_fit;
    lo = std::min<double>(lo, e.self_fit);
    hi = std::max<double>(hi, e.self_fit);
  }
  std::printf(
    "self-fit written for %zu/%zu entries: min %.3f mean %.3f max %.3f\n", done, db.entries.size(),
    lo, sum / static_cast<double>(db.entries.size()), hi);
  return 0;
}
