#pragma once

#include <vector>
#include <tuple>
#include <small_gicp/util/lie.hpp>
#include <small_gicp/factors/linearization.hpp>

namespace small_gicp {

/// @brief Single-thread reduction
struct SerialReduction {
  /// @brief Sum up linearized systems
  /// @param target       Target point cloud
  /// @param source       Source point cloud
  /// @param target_tree  Nearest neighbor search for the target point cloud
  /// @param rejector     Correspondence rejector
  /// @param T            Linearization point
  /// @param factors      Factors to be linearized
  /// @return             Sum of the linearized systems (information matrix, information vector, and error)
  template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree, typename CorrespondenceRejector, typename Factor>
  std::tuple<factor_detail::Hessian, factor_detail::Gradient, double> linearize(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const CorrespondenceRejector& rejector,
    const Transform& T,
    std::vector<Factor>& factors) const {
    using earth::mantle::collection::Matrix;
    using earth::mantle::math::operator+;
    Matrix<double, 6, 6> sum_H{};
    Matrix<double, 6, 1> sum_b{};
    double sum_e = 0.0;

    for (size_t i = 0; i < factors.size(); i++) {
      factor_detail::Hessian H{};
      factor_detail::Gradient b{};
      double e;

      if (!factors[i].linearize(target, source, target_tree, T, i, rejector, &H, &b, &e)) {
        continue;
      }

      sum_H = sum_H + H;
      sum_b = sum_b + b;
      sum_e += e;
    }

    return {sum_H, sum_b, sum_e};
  }

  /// @brief Sum up evaluated errors
  /// @param target       Target point cloud
  /// @param source       Source point cloud
  /// @param T            Evaluation point
  /// @param factors      Factors to be evaluated
  /// @return Sum of the evaluated errors
  template <typename TargetPointCloud, typename SourcePointCloud, typename Factor>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T, std::vector<Factor>& factors) const {
    double sum_e = 0.0;
    for (size_t i = 0; i < factors.size(); i++) {
      sum_e += factors[i].error(target, source, T);
    }
    return sum_e;
  }
};

}  // namespace small_gicp
