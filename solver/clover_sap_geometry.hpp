#pragma once
#include <stdexcept>
#include <vector>

namespace clover_qxs {
// The packed coarse lattice needs a full SIMD vector in x and y. Use
// four-site SAP blocks when possible and the two-site kernel otherwise.
inline std::vector<int> sap_blocks(const std::vector<int>& local, int vx, int vy) {
  if (local.size() != 4) throw std::invalid_argument("QXS requires four local extents");
  const int lanes[] = {vx, vy, 1, 1};
  std::vector<int> blocks(4, 4);
  for (int axis = 0; axis < 4; ++axis) {
    if (axis < 2 && local[axis] % (4 * lanes[axis]) != 0) blocks[axis] = 2;
    if (local[axis] <= 0 || local[axis] % (blocks[axis] * lanes[axis]) != 0)
      throw std::invalid_argument("local lattice is incompatible with QXS SAP and coarse packing");
  }
  return blocks;
}
}  // namespace clover_qxs
