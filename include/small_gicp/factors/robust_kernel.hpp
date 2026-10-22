// SPDX-FileCopyrightText: Copyright 2024 Kenji Koide
// SPDX-License-Identifier: MIT
#pragma once

#include <small_gicp/util/lie.hpp>
#include <small_gicp/factors/linearization.hpp>

namespace small_gicp {

/// @brief Huber robust kernel
struct Huber {
public:
  /// @brief Huber robust kernel setting
  struct Setting {
    double c = 1.0;  ///< Kernel width
  };

  /// @brief Constructor
  Huber(const Setting& setting) : c(setting.c) {}

  /// @brief Compute the weight for an error
  /// @param e Error
  /// @return  Weight
  double weight(double e) const {
    const double e_abs = std::abs(e);
    return e_abs < c ? 1.0 : c / e_abs;
  }

public:
  const double c;  ///< Kernel width
};

/// @brief Cauchy robust kernel
struct Cauchy {
public:
  /// @brief Huber robust kernel setting
  struct Setting {
    double c = 1.0;  ///< Kernel width
  };

  /// @brief Constructor
  Cauchy(const Setting& setting) : c(setting.c) {}

  /// @brief Compute the weight for an error
  /// @param e Error
  /// @return  Weight
  double weight(double e) const { return c / (c + e * e); }

public:
  const double c;  ///< Kernel width
};

/// @brief Robustify a factor with a robust kernel
/// @tparam Kernel  Robust kernel
/// @tparam Factor  Factor
template <typename Kernel, typename Factor>
struct RobustFactor {
public:
  /// @brief Robust factor setting
  struct Setting {
    typename Kernel::Setting robust_kernel;  ///< Robust kernel setting
    typename Factor::Setting factor;         ///< Factor setting
  };

  /// @brief Constructor
  RobustFactor(const Setting& setting = Setting()) : robust_kernel(setting.robust_kernel), factor(setting.factor) {}

  /// @brief Linearize the factor
  template <typename TargetPointCloud, typename SourcePointCloud, typename TargetTree, typename CorrespondenceRejector>
  bool linearize(
    const TargetPointCloud& target,
    const SourcePointCloud& source,
    const TargetTree& target_tree,
    const Transform& T,
    size_t source_index,
    const CorrespondenceRejector& rejector,
    factor_detail::Hessian* H,
    factor_detail::Gradient* b,
    double* e) {
    if (!factor.linearize(target, source, target_tree, T, source_index, rejector, H, b, e)) {
      return false;
    }

    // Robustify the linearized factor
    const double w = robust_kernel.weight(std::sqrt(*e));
    using earth::mantle::math::operator*;
    *H = *H * w;
    *b = *b * w;
    *e *= w;

    return true;
  }

  /// @brief Evaluate error
  template <typename TargetPointCloud, typename SourcePointCloud>
  double error(const TargetPointCloud& target, const SourcePointCloud& source, const Transform& T) const {
    const double e = factor.error(target, source, T);
    return robust_kernel.weight(std::sqrt(e)) * e;
  }

  /// @brief  Check if the factor is inlier
  bool inlier() const { return factor.inlier(); }

public:
  Kernel robust_kernel;  ///< Robust kernel
  Factor factor;         ///< Factor
};

}  // namespace small_gicp
