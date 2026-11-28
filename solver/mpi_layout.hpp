#pragma once
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace clover_layout {
inline bool compatible(const std::vector<int>& shape,const std::vector<int>& grid,
                       const std::vector<int>& packing) {
  for(int d=0;d<4;++d)
    if(shape[d]<=0 || grid[d]<=0 || shape[d]%grid[d] || (shape[d]/grid[d])%packing[d])
      return false;
  return true;
}
inline std::vector<int> divisors(int value) {
  std::vector<int> result;
  for(int k=1;k<=value/k;++k) if(value%k==0) {
    result.push_back(k);
    if(k!=value/k) result.push_back(value/k);
  }
  std::sort(result.begin(),result.end());
  return result;
}

// For equal local volumes, summed face/volume is proportional to halo traffic.
// Ties keep x/y as contiguous as possible for the packed kernels.
inline std::vector<int> choose(const std::vector<int>& shape,int ranks,const std::vector<int>& packing) {
  if(shape.size()!=4 || packing.size()!=4 || ranks<=0)
    throw std::invalid_argument("invalid MPI layout request");
  std::vector<int> best;
  long double lowest=std::numeric_limits<long double>::infinity();
  for(int x:divisors(ranks)) for(int y:divisors(ranks/x)) for(int z:divisors(ranks/x/y)) {
    const std::vector<int> grid{x,y,z,ranks/x/y/z};
    long double cost=0;
    if(!compatible(shape,grid,packing)) continue;
    for(int d=0;d<4;++d) {
      if(grid[d]>1) cost+=static_cast<long double>(grid[d])/shape[d];
    }
    if(cost<lowest) {best=grid;lowest=cost;}
  }
  return best;
}
}
