// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <tbb/tbb.h>
#include <small_gicp/util/lie.hpp>
#include <tuple>
#include <vector>
#include <small_gicp/factors/linearization.hpp>

namespace small_gicp {

/// @brief Summation for linearized systems
template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree, typename CorrespondenceRejector, typename Factor>
struct LinearizeSum {
  LinearizeSum(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const CorrespondenceRejector& rejector,
    const Transform& T,
    std::vector<Factor>& factors)
  : target(target),
    source(source),
    target_tree(target_tree),
    rejector(rejector),
    T(T),
    factors(factors),
    H{},
    b{},
    e(0.0) {}

  LinearizeSum(LinearizeSum& x, tbb::split)
  : target(x.target),
    source(x.source),
    target_tree(x.target_tree),
    rejector(x.rejector),
    T(x.T),
    factors(x.factors),
    H{},
    b{},
    e(0.0) {}

  void operator()(const tbb::blocked_range<size_t>& r) {
    using earth::mantle::math::operator+;
    factor_detail::Hessian Ht = H;
    factor_detail::Gradient bt = b;
    double et = e;

    for (size_t i = r.begin(); i != r.end(); i++) {
      factor_detail::Hessian Hi;
      factor_detail::Gradient bi;
      double ei;

      if (!factors[i].linearize(target, source, target_tree, T, i, rejector, &Hi, &bi, &ei)) {
        continue;
      }

      Ht = Ht + Hi;
      bt = bt + bi;
      et += ei;
    }

    H = Ht;
    b = bt;
    e = et;
  }

  void join(const LinearizeSum& y) {
    using earth::mantle::math::operator+;
    H = H + y.H;
    b = b + y.b;
    e += y.e;
  }

  const TargetPointCloud& target;
  const SourcePointCloud& source;
  const TargetTree& target_tree;
  const CorrespondenceRejector& rejector;
  const Transform& T;
  std::vector<Factor>& factors;

  factor_detail::Hessian H;
  factor_detail::Gradient b;
  double e;
};

/// @brief Summation for evaluated errors
template <typename TargetPointCloud, typename SourcePointCloud, typename Factor>
struct ErrorSum {
  ErrorSum(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T, std::vector<Factor>& factors)
  : target(target),
    source(source),
    T(T),
    factors(factors),
    e(0.0) {}

  ErrorSum(ErrorSum& x, tbb::split) : target(x.target), source(x.source), T(x.T), factors(x.factors), e(0.0) {}

  void operator()(const tbb::blocked_range<size_t>& r) {
    double et = e;
    for (size_t i = r.begin(); i != r.end(); i++) {
      et += factors[i].error(target, source, T);
    }
    e = et;
  }

  void join(const ErrorSum& y) { e += y.e; }

  const TargetPointCloud& target;
  const SourcePointCloud& source;
  const Transform& T;
  std::vector<Factor>& factors;

  double e;
};

/// @brief Parallel reduction with TBB backend
struct ParallelReductionTBB {
  ParallelReductionTBB() {}

  template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree, typename CorrespondenceRejector, typename Factor>
  std::tuple<factor_detail::Hessian, factor_detail::Gradient, double> linearize(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const CorrespondenceRejector& rejector,
    const Transform& T,
    std::vector<Factor>& factors) const {
    //
    LinearizeSum<TargetPointCloud, SourcePointCloud, TargetTree, CorrespondenceRejector, Factor> sum(target, source, target_tree, rejector, T, factors);

    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, factors.size(), 8), sum);

    return {sum.H, sum.b, sum.e};
  }

  template <typename TargetPointCloud, typename SourcePointCloud, typename Factor>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T, std::vector<Factor>& factors) const {
    ErrorSum<TargetPointCloud, SourcePointCloud, Factor> sum(target, source, T, factors);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, factors.size(), 16), sum);
    return sum.e;
  }
};

}  // namespace small_gicp
