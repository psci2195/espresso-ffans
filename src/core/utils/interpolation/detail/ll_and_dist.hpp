#ifndef UTILS_INTERPOLATION_DETAIL_LL_AND_DIST_HPP
#define UTILS_INTERPOLATION_DETAIL_LL_AND_DIST_HPP

#include "Vector.hpp"

#include "pos_shift.hpp"

#include <array>
#include <cmath>
#include <utility>

namespace Utils {
namespace Interpolation {
namespace detail {
template <size_t order>
std::pair<std::array<int, 3>, std::array<double, 3>>
ll_and_dist(const Vector3d &pos, const Vector3d &grid_spacing,
            const Vector3d &offset) {
  /* Distance to the nearest mesh point in units of h \in [-0.5, 0.5) */
  std::array<double, 3> dist;
  /* Index of the lower left corner of the assignment cube */
  std::array<int, 3> ll;

  for (int dim = 0; dim < 3; dim++) {
    const double nmp_pos = (pos[dim] - offset[dim]) / grid_spacing[dim] +
                           detail::pos_shift<order>();
    const double nmp_ind = static_cast<int>(nmp_pos);
    dist[dim] = nmp_pos - nmp_ind - 0.5;
    ll[dim] = nmp_ind - (order - 1) / 2;
  }

  return {ll, dist};
}
}
}
}

#endif
