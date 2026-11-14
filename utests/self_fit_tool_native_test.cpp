#include <alg/point_cloud_io.h>
#include <alg/scan_context.h>
#include <robot_core/test/doctest.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace scan_context = earth::mantle::alg::scan_context;

namespace
{
using Cloud = earth::mantle::alg::PointCloud<earth::mantle::alg::PointXYZI>;

struct Fixture
{
  std::string directory;
  std::string map, database, scans;
  Fixture()
  {
    char pattern[] = "/tmp/sc_self_fit_test.XXXXXX";
    const char * path = mkdtemp(pattern);
    if (!path) return;
    directory = path;
    map = directory + "/map.pcd";
    database = directory + "/map.sc";
    scans = directory + "/map.scans";
  }
  ~Fixture()
  {
    if (directory.empty()) return;
    unlink(map.c_str());
    unlink(database.c_str());
    unlink(scans.c_str());
    rmdir(directory.c_str());
  }
};

int runTool(const Fixture & fixture, const char * leaf = nullptr)
{
  const pid_t child = fork();
  if (child == -1) return -1;
  if (child == 0) {
    if (leaf)
      setenv("SF_LEAF", leaf, 1);
    else
      unsetenv("SF_LEAF");
    execl(
      SELF_FIT_TOOL_PATH, SELF_FIT_TOOL_PATH, fixture.map.c_str(), fixture.database.c_str(),
      fixture.scans.c_str(), static_cast<char *>(nullptr));
    _exit(127);
  }
  int status = 0;
  if (waitpid(child, &status, 0) != child || !WIFEXITED(status)) return -1;
  return WEXITSTATUS(status);
}

std::string contents(const std::string & path)
{
  std::ifstream file(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}
}  // namespace

TEST_CASE("self-fit tool scores matching mismatching and unusable clouds")
{
  Fixture fixture;
  REQUIRE_FALSE(fixture.directory.empty());
  Cloud map;
  std::vector<std::vector<float>> scans(6);
  for (int y = 0; y <= 30; ++y)
    for (int z = 0; z <= 20; ++z) {
      const float py = float(y) * .2f, pz = float(z) * .2f;
      map.push_back({0, py, pz, 0});
      scans[0].insert(scans[0].end(), {0, py, pz});
      scans[1].insert(scans[1].end(), {2, py, pz});
      scans[2].insert(scans[2].end(), {y < 15 ? 0.f : 2.f, py, pz});
      scans[3].insert(scans[3].end(), {py, pz, 0});  // horizontal normal: excluded
    }
  scans[4] = {0, 0, 0};  // below minimum sampled cloud size
  REQUIRE(earth::mantle::alg::writePCDBinary(fixture.map, map) == 0);
  scan_context::Database database;
  for (std::size_t i = 0; i < scans.size(); ++i)
    REQUIRE(database.add(
      0, 0, 0, 0, std::vector<float>(database.cfg.num_ring * database.cfg.num_sector, 1), 0));
  REQUIRE(database.save(fixture.database));
  REQUIRE(scan_context::saveScans(fixture.scans, scans));
  REQUIRE(runTool(fixture) == 0);
  scan_context::Database scored;
  REQUIRE(scored.load(fixture.database));
  REQUIRE(scored.entries.size() == scans.size());
  CHECK(scored.entries[0].self_fit == doctest::Approx(1));
  CHECK(scored.entries[1].self_fit == 0);
  CHECK(scored.entries[2].self_fit > 0);
  CHECK(scored.entries[2].self_fit < 1);
  for (std::size_t i : {3u, 4u, 5u}) CHECK(scored.entries[i].self_fit == 0);
  for (const auto & entry : scored.entries) {
    CHECK(entry.x == 0);
    CHECK(entry.y == 0);
    CHECK(entry.z == 0);
    CHECK(entry.yaw == 0);
    CHECK(entry.desc == database.entries[0].desc);
  }
  const auto before = contents(fixture.database);
  for (const char * leaf : {"0", "-1", "nan", "inf"}) {
    CHECK(runTool(fixture, leaf) == 1);
    CHECK(contents(fixture.database) == before);
  }
  for (const std::string & truncated :
       {std::string("SCSC"), std::string("SCSC0001"), std::string("SCSC0001\1", 9)}) {
    {
      std::ofstream file(fixture.scans, std::ios::binary);
      file.write(truncated.data(), truncated.size());
    }
    CHECK(runTool(fixture) == 1);
    CHECK(contents(fixture.database) == before);
  }
}
