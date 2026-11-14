#include <random>
#include <algorithm>
#include <small_gicp/util/sort_omp.hpp>
#include <small_gicp/util/sort_tbb.hpp>

#include <robot_core/test/doctest.h>
#include <cmath>

using namespace small_gicp;

// Check if two vectors are identical
template <typename T>
bool identical(const std::vector<T>& arr1, const std::vector<T>& arr2) {
  if (arr1.size() != arr2.size()) {
    return false;
  }

  for (size_t i = 0; i < arr1.size(); i++) {
    if (arr1[i] != arr2[i]) {
      return false;
    }
  }
  return true;
}

template <typename T>
void test_radix_sort(std::mt19937& mt) {
  std::uniform_int_distribution<> size_dist(0, 8192);

  for (int i = 0; i < 20; i++) {
    std::vector<T> data(size_dist(mt));
    std::generate(data.begin(), data.end(), [&] { return mt(); });

    std::vector<T> sorted = data;
    std::stable_sort(sorted.begin(), sorted.end());

    std::vector<T> sorted_tbb = data;
    radix_sort_tbb(sorted_tbb.data(), sorted_tbb.data() + sorted_tbb.size(), [](const T x) { return x; });

    INFO("i=", i, ", N=", data.size());
    CHECK(identical(sorted, sorted_tbb));
  }

  for (int i = 0; i < 20; i++) {
    std::vector<std::pair<T, std::uint64_t>> data(size_dist(mt));
    std::generate(data.begin(), data.end(), [&] { return std::make_pair<T, std::uint64_t>(mt(), mt()); });

    std::vector<std::pair<T, std::uint64_t>> sorted = data;
    std::stable_sort(sorted.begin(), sorted.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

    std::vector<std::pair<T, std::uint64_t>> sorted_tbb = data;
    radix_sort_tbb(sorted_tbb.data(), sorted_tbb.data() + sorted_tbb.size(), [](const auto& x) -> T { return x.first; });

    INFO("i=", i, ", N=", data.size());
    CHECK(identical(sorted, sorted_tbb));
  }
}

// Test radix_sort_tbb
TEST_CASE("SortTBB.RadixSortTest") {
  std::mt19937 mt;

  test_radix_sort<std::uint8_t>(mt);
  test_radix_sort<std::uint16_t>(mt);
  test_radix_sort<std::uint32_t>(mt);
  test_radix_sort<std::uint64_t>(mt);

  // empty
  std::vector<std::uint64_t> empty_vector;
  radix_sort_tbb(empty_vector.data(), empty_vector.data() + empty_vector.size(), [](const std::uint64_t x) { return x; });
  CHECK_MESSAGE((empty_vector.empty()), ("Empty vector check"));
}
