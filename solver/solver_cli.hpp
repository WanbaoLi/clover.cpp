#pragma once
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct SolveOptions {
  bool plan_layout=false, checkpoint=false, compact_eo=false, profile_mg=false;
  std::string resume_from, precision="double", mg_outer="fbicgstab";
  std::string gauge, output, source = "point", fields = "full", method = "auto";
  std::vector<int> lattice, grid, position{0,0,0,0}, rhs, sap_block;
  double mass = -0.279, csw = 1.160678675188, rho = 0.125;
  double tolerance = 5e-13, acceptance = 1e-12;
  int bicg_iterations = 0, bicg_cycles = 0;  // Zero selects the method-specific default.
  int mg_max_iters = 0, gcr_restart = 16;
  int threads = 1, vectors = 24, setup_steps = 4, smooth_steps = 4;
};

inline const char* solve_usage() {
  return "clover-solve --gauge FILE|unit --lattice X,Y,Z,T [--grid auto|PX,PY,PZ,PT] --output NEW_DIRECTORY\n"
         "  clover-solve --plan-layout --lattice X,Y,Z,T [--grid auto|PX,PY,PZ,PT] [--method auto|mg|bicgstab|eo-bicgstab]\n"
         "  [--method auto|mg|bicgstab|eo-bicgstab] [--mass M] [--csw C] [--rho R] [--source point|wall] [--position X,Y,Z,T]\n"
         "  [--rhs all|0,1,...,11] [--fields full|none] [--tol R] [--accept R]\n"
         "  [--bicg-iters N] [--bicg-cycles N] [--compact-eo] [--checkpoint] [--resume-from DIRECTORY]\n"
         "  [--threads N] [--nvec N] [--setup-steps N] [--smooth-iters N] [--profile-mg]\n"
         "  [--sap-block X,Y,Z,T] [--mg-max-iters N] [--precision single|double]\n"
         "  [--mg-outer fbicgstab|gcr] [--gcr-restart N]\n";
}
inline int solve_int(const std::string& s) {
  char* end=nullptr;
  errno=0;
  const long v=std::strtol(s.c_str(),&end,10);
  if (s.empty() || *end || errno || v<std::numeric_limits<int>::min() || v>std::numeric_limits<int>::max())
    throw std::runtime_error("invalid integer: "+s);
  return static_cast<int>(v);
}
inline double solve_double(const std::string& s) {
  char* end=nullptr;
  errno=0;
  const double v=std::strtod(s.c_str(),&end);
  if (s.empty() || *end || errno || !std::isfinite(v))
    throw std::runtime_error("invalid finite number: "+s);
  return v;
}
inline std::vector<int> solve_list(const std::string& s) {
  if (s.empty() || s.back()==',') throw std::runtime_error("empty list item: "+s);
  std::istringstream input(s);
  std::string item;
  std::vector<int> result;
  while (std::getline(input,item,',')) result.push_back(solve_int(item));
  return result;
}
inline SolveOptions solve_options(int argc,char** argv) {
  SolveOptions o;
  const char* env=std::getenv("OMP_NUM_THREADS");
  if (env) o.threads=solve_int(env);
  std::set<std::string> seen;
  for (int i=1;i<argc;++i) {
    const std::string k=argv[i];
    if (!seen.insert(k).second) throw std::runtime_error("duplicate option: "+k);
    if (k=="--plan-layout") {o.plan_layout=true;continue;}
    if (k=="--profile-mg") {o.profile_mg=true;continue;}
    if (k=="--compact-eo") {o.compact_eo=true;continue;}
    if (k=="--checkpoint") {o.checkpoint=true;continue;}
    if (++i>=argc) throw std::runtime_error("missing value: "+k);
    const std::string v=argv[i];
    if (k=="--gauge") o.gauge=v;
    else if (k=="--output") o.output=v;
    else if (k=="--resume-from") {
      if(v.empty()) throw std::runtime_error("--resume-from requires a directory");
      o.resume_from=v;o.checkpoint=true;
    }
    else if (k=="--lattice") o.lattice=solve_list(v);
    else if (k=="--grid") { if(v!="auto") o.grid=solve_list(v); }
    else if (k=="--position") o.position=solve_list(v);
    else if (k=="--source") o.source=v;
    else if (k=="--precision") {
      if (v!="single" && v!="double") throw std::runtime_error("--precision must be single or double");
      o.precision=v;
    }
    else if (k=="--method") o.method=v;
    else if (k=="--fields") o.fields=v;
    else if (k=="--mass") o.mass=solve_double(v);
    else if (k=="--csw") o.csw=solve_double(v);
    else if (k=="--rho") o.rho=solve_double(v);
    else if (k=="--tol") o.tolerance=solve_double(v);
    else if (k=="--accept") o.acceptance=solve_double(v);
    else if (k=="--bicg-iters" || k=="--bicg-cycles") {
      const int n=solve_int(v);
      if (n<1) throw std::runtime_error(k+" must be positive");
      if (k=="--bicg-iters") o.bicg_iterations=n;
      else o.bicg_cycles=n;
    }
    else if (k=="--threads") o.threads=solve_int(v);
    else if (k=="--sap-block") {
      o.sap_block=solve_list(v);
      if (o.sap_block.size()!=4 ||
          std::any_of(o.sap_block.begin(),o.sap_block.end(),[](int n) { return n<2 || n%2; }))
        throw std::runtime_error("--sap-block requires four positive even extents");
    }
    else if (k=="--mg-max-iters") {
      o.mg_max_iters=solve_int(v);
      if (o.mg_max_iters<1 || o.mg_max_iters>(std::numeric_limits<int>::max()-1)/2)
        throw std::runtime_error("--mg-max-iters out of range");
    }
    else if (k=="--mg-outer") {
      if (v!="fbicgstab" && v!="gcr") throw std::runtime_error("--mg-outer must be fbicgstab or gcr");
      o.mg_outer=v;
    }
    else if (k=="--gcr-restart") {
      o.gcr_restart=solve_int(v);
      if (o.gcr_restart<1) throw std::runtime_error("--gcr-restart must be positive");
    }
    else if (k=="--nvec") o.vectors=solve_int(v);
    else if (k=="--setup-steps") o.setup_steps=solve_int(v);
    else if (k=="--smooth-iters") o.smooth_steps=solve_int(v);
    else if (k=="--rhs") { if (v!="all") o.rhs=solve_list(v); }
    else throw std::runtime_error("unknown option: "+k);
  }
  if (o.rhs.empty()) for (int c=0;c<12;++c) o.rhs.push_back(c);
  if ((!o.plan_layout && (o.gauge.empty() || o.output.empty())) || o.lattice.size()!=4 ||
      (!o.grid.empty() && o.grid.size()!=4) || o.position.size()!=4)
    throw std::runtime_error("supply four lattice extents and a four-dimensional or automatic grid; solving requires gauge and output");
  if ((o.fields!="full" && o.fields!="none") || (o.source!="point" && o.source!="wall"))
    throw std::runtime_error("fields must be full/none and source must be point/wall");
  if (o.method!="auto" && o.method!="mg" && o.method!="bicgstab" && o.method!="eo-bicgstab") throw std::runtime_error("method must be auto, mg, bicgstab or eo-bicgstab");
  if (o.mass<=-4 || o.tolerance<=0 || o.tolerance>=1 || o.tolerance*o.tolerance==0 || o.acceptance<=0 || o.acceptance>=1 ||
      o.threads<1 || o.vectors<1 || o.setup_steps<1 || o.smooth_steps<1)
    throw std::runtime_error("invalid mass, tolerance, thread count, or multigrid parameters");
  if (std::set<int>(o.rhs.begin(),o.rhs.end()).size()!=o.rhs.size() || *std::min_element(o.rhs.begin(),o.rhs.end())<0 || *std::max_element(o.rhs.begin(),o.rhs.end())>=12)
    throw std::runtime_error("source columns must be distinct integers from 0 through 11");
  if (seen.count("--gcr-restart") && o.mg_outer!="gcr")
    throw std::runtime_error("--gcr-restart requires --mg-outer gcr");
  if (o.compact_eo && o.method!="eo-bicgstab") throw std::runtime_error("--compact-eo requires --method eo-bicgstab");
  return o;
}
inline std::string solve_join(const std::vector<int>& v) {
  std::ostringstream s;
  for (size_t i=0;i<v.size();++i) { if (i) s<<','; s<<v[i]; }
  return s.str();
}
