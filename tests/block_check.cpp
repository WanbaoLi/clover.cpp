#include "bridge_setup.h"
#include "bridge_init_factory.h"
#include "lib/IO/gaugeConfig.h"
#include "lib/Fopr/fopr.h"
#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <mpi.h>

// The scalar Wilson difference replaces only hopping links at block faces.
// Clover's diagonal is evaluated on all original links.
int main(int argc,char** argv) {
  if(argc!=1 && argc!=3) return 2;
  const bool from_file=argc==3;
  std::string file=from_file ? argv[1] : "";
  int gt=from_file ? std::atoi(argv[2]) : 0;
  bridge_initialize(&argc,&argv);
  int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  if(!from_file) gt=4*ranks;
  int result=1;
  {
    const int threads=std::getenv("OMP_NUM_THREADS") ? std::atoi(std::getenv("OMP_NUM_THREADS")) : 1;
    const int xy=from_file ? 24 : 8;
    bridge_setup({xy,xy,4,gt},{1,1,1,ranks},threads,3,"stdout","Crucial");
    const int n=CommonParameters::Nvol();
    const std::vector<int> block{2,2,4,4};
    const int length[]={xy,xy,4,gt/ranks};
    Field_G u(n,4),cut(n,4);
    if(from_file) GaugeConfig("ILDG_Parallel").read(u,file);
    else {
      // Product of two noncommuting SU(3) rotations. The global site index
      // fixes the same input independently of MPI rank numbering.
      for(int s=0;s<n;++s) for(int d=0;d<4;++d) {
        int global_site=s+Communicator::ipe(3)*n;
        double a=.21*std::sin(.13*global_site+.28*d);
        double b=.18*std::cos(.09*global_site+.31*d);
        double ca=std::cos(a),sa=std::sin(a),cb=std::cos(b),sb=std::sin(b);
        std::complex<double> row[]={{ca,0},{0,sa*cb},{0,sa*sb},
                                    {0,sa},{ca*cb,0},{ca*sb,0},
                                    {0,0},{-sb,0},{cb,0}};
        for(int k=0;k<9;++k) u.set_ri(k,s,d,row[k].real(),row[k].imag());
      }
    }
    copy(cut,u);
    Field_F b,ref,full,tmp,check,expected;
    for(int s=0;s<n;++s) {
      int site=s;
      for(int d=0;d<4;++d) {
        const int pos=site%length[d];site/=length[d];
        if((pos+1)%block[d]==0)
          for(int c=0;c<9;++c) cut.set_ri(c,s,d,0.,0.);
      }
      for(int j=0;j<24;++j) b.ptr(24*s)[j]=std::sin(.09*double(j+17*s+3*rank));
    }
    Parameters p;
    p.set_string("fermion_type","Clover");p.set_string("gamma_matrix_type","Dirac");
    p.set_double("hopping_parameter",.5/(4.-.279));p.set_double("clover_coefficient",1.160678675188);
    p.set_int_vector("boundary_condition",{1,1,1,-1});p.set_int_vector("block_size",block);
    p.set_string("verbose_level","Crucial");
    std::unique_ptr<Fopr> c(Fopr::New("Clover",p)),w(Fopr::New("Wilson",p));
    c->set_mode("D");w->set_mode("D");c->set_config(&u);w->set_config(&u);
    c->mult(full,b);copy(ref,full);w->mult(tmp,b);axpy(ref,-1.,tmp);
    w->set_config(&cut);w->mult(tmp,b);axpy(ref,1.,tmp);
    using Packed=AField<float,QXS>;
    AFopr_Clover_dd<Packed> dd(p);dd.set_config(&u);dd.set_mode("D");
    Packed x(24,n,1),out(24,n,1),hop(24,n,1),sum(24,n,1);
    double max_error=0;
    auto measure=[&](const char* name,const Field_F& want) {
      axpy(check,-1.,want);
      const double err=std::sqrt(check.norm2()/want.norm2());
      if(!rank) std::printf("BLOCK_CHECK %s %.17g\n",name,err);
      if(!std::isfinite(err)) max_error=1;
      else max_error=std::max(max_error,err);
    };
#pragma omp parallel
    {dd.convert(x,b);dd.mult_dd(sum,x);dd.reverse(check,sum);}
    measure("Dirichlet",ref);
    for(int eo=0;eo<2;++eo) {
      copy(expected,ref);
      for(int site=0;site<n;++site) {
        int v=site;
        int parity=Communicator::ipe(3)*(length[3]/block[3]);
        for(int d=0;d<4;++d){parity+=(v%length[d])/block[d];v/=length[d];}
        if(parity%2!=eo) std::fill_n(expected.ptr(24*site),24,0.);
      }
#pragma omp parallel
      {out.set(0.);dd.mult_sap(out,x,eo);dd.reverse(check,out);}
      measure(eo?"SAP_odd":"SAP_even",expected);
    }
#pragma omp parallel
    {
      for(int d=0;d<4;++d) {
        dd.mult_dup(hop,x,d);axpy(sum,1.0f,hop);
        dd.mult_ddn(hop,x,d);axpy(sum,1.0f,hop);
      }
      dd.reverse(check,sum);
    }
    measure("reconstructed_full_D",full);
    result=max_error<2e-6 ? 0 : 1;
  }
  bridge_finalize();return result;
}
