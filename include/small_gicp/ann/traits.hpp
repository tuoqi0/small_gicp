#pragma once

#include <type_traits>
#include <cstddef>
#include <utility>
#include <small_gicp/points/traits.hpp>

namespace small_gicp {

namespace ann_detail {
template <typename Point, typename Query>
double squared_distance(const Point& point, const Query& query) {
  const double dx = traits::detail::coordinate(point, 0) - traits::detail::coordinate(query, 0);
  const double dy = traits::detail::coordinate(point, 1) - traits::detail::coordinate(query, 1);
  const double dz = traits::detail::coordinate(point, 2) - traits::detail::coordinate(query, 2);
  const double dw = traits::detail::coordinate(point, 3) - traits::detail::coordinate(query, 3);
  // Preserve the former four-coordinate packet reduction order.
  return (dx * dx + dz * dz) + (dy * dy + dw * dw);
}
}

namespace traits {

template <typename T>
struct Traits;

/// @brief Find k-nearest neighbors.
/// @param tree       Nearest neighbor search (e.g., KdTree)
/// @param point      Query point
/// @param k          Number of neighbors
/// @param k_indices  [out] Indices of k-nearest neighbors
/// @param k_sq_dists [out] Squared distances to k-nearest neighbors
/// @return Number of found neighbors
template <typename T, typename Point>
size_t knn_search(const T& tree, const Point& point, size_t k, size_t* k_indices, double* k_sq_dists) {
  return Traits<T>::knn_search(tree, point, k, k_indices, k_sq_dists);
}

/// @brief Check if T has nearest_neighbor_search method.
template <typename T, typename Point, typename = void>
struct has_nearest_neighbor_search : std::false_type {};

template <typename T, typename Point>
struct has_nearest_neighbor_search<T, Point, std::void_t<decltype(
  Traits<T>::nearest_neighbor_search(std::declval<const T&>(),
    std::declval<const Point&>(), std::declval<size_t*>(), std::declval<double*>()))>>
  : std::true_type {};

/// @brief Find the nearest neighbor. If Traits<T>::nearest_neighbor_search is not defined, fallback to knn_search with k=1.
/// @param tree       Nearest neighbor search (e.g., KdTree)
/// @param point      Query point
/// @param k_index    [out] Index of the nearest neighbor
/// @param k_sq_dist  [out] Squared distance to the nearest neighbor
/// @return 1 if a neighbor is found else 0
template <typename T, typename Point,
          std::enable_if_t<has_nearest_neighbor_search<T, Point>::value, bool> = true>
size_t nearest_neighbor_search(const T& tree, const Point& point, size_t* k_index, double* k_sq_dist) {
  return Traits<T>::nearest_neighbor_search(tree, point, k_index, k_sq_dist);
}

/// @brief Find the nearest neighbor. If Traits<T>::nearest_neighbor_search is not defined, fallback to knn_search with k=1.
/// @param tree       Nearest neighbor search (e.g., KdTree)
/// @param point      Query point
/// @param k_index    [out] Index of the nearest neighbor
/// @param k_sq_dist  [out] Squared distance to the nearest neighbor
/// @return 1 if a neighbor is found else 0
template <typename T, typename Point,
          std::enable_if_t<!has_nearest_neighbor_search<T, Point>::value, bool> = true>
size_t nearest_neighbor_search(const T& tree, const Point& point, size_t* k_index, double* k_sq_dist) {
  return Traits<T>::knn_search(tree, point, 1, k_index, k_sq_dist);
}

}  // namespace traits

}  // namespace small_gicp
