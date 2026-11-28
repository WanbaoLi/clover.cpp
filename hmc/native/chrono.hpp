#pragma once
#include "schur.hpp"
#include <deque>
#include <map>
#include <utility>

// Least-residual chronological predictor in the real inner product on complex fields.
// Recompute each image with the current gauge. Distinct mass pairs keep separate histories.
using SolutionHistory=std::deque<Field_F>;
static std::map<std::pair<double,double>,SolutionHistory> past_solutions;

inline int chronological_guess(OddSchur& M,const SolutionHistory& stored,const Field_F& rhs,
                               Field_F& x,Field_F& residual) {
  std::vector<Field_F> directions,images;
  Field_F y(half_volume,1);
  const double rhs_norm=rhs.norm();
  for(auto prev=stored.rbegin();prev!=stored.rend();++prev) {
    Field_F v(*prev),a(half_volume,1);
    M.apply(y,v,false);M.apply(a,y,true);
    const double before=a.norm();
    if(!(before>0.)) continue;
    // Modified Gram-Schmidt twice for nearly collinear neighboring MD solutions.
    for(int pass=0;pass<2;++pass) for(size_t j=0;j<images.size();++j) {
      const double coefficient=dot(images[j],a);
      axpy(a,-coefficient,images[j]);axpy(v,-coefficient,directions[j]);
    }
    const double length=a.norm();
    if(length<1.e-10*before) continue;
    scal(a,1./length);scal(v,1./length);
    const double coefficient=dot(a,rhs);
    axpy(x,coefficient,v);
    images.push_back(a);directions.push_back(v);
  }
  // Evaluate b-MdagM*x explicitly before CG; a weak or degenerate prediction uses zero.
  M.apply(y,x,false);M.apply(residual,y,true);aypx(-1.,residual,rhs);
  const double predicted=residual.norm();
  if(!std::isfinite(predicted) || predicted>rhs_norm) {
    x.set(0.);copy(residual,rhs);return 0;
  }
  return static_cast<int>(images.size());
}
