#if defined(USE_ALT_QXS) && defined(USE_GROUP_SU3) && defined(__ARM_FEATURE_SVE_BITS) && __ARM_FEATURE_SVE_BITS == 512
#include <arm_sve.h>
#include <omp.h>
#include <stdexcept>

// FP64 Bridge++ QXS layout: each component is one vector of eight different
// lattice sites. The leading RHS dimension stays in separate existing fields.
namespace {
using Vec = svfloat64_t __attribute__((arm_sve_vector_bits(512)));
constexpr int lanes = sizeof(Vec)/sizeof(double);
static_assert(lanes == 8, "This kernel requires the FP64 512-bit QXS packing");
constexpr int spin_components = 24;
constexpr int clover_components = 72;
// Upper triangle of the six-by-six Hermitian block in set_csw_chrot().
constexpr int packed[6][6] = {
    {0,6,16,28,24,34}, {-1,1,12,18,32,26}, {-1,-1,2,8,20,30},
    {-1,-1,-1,3,14,22}, {-1,-1,-1,-1,4,10}, {-1,-1,-1,-1,-1,5}
};

// Every matrix coefficient load precedes the unrolled RHS loop. Each row's
// four complex accumulators remain live across the six columns.
template<int Row, int Batch>
inline void matvec_row(const double* matrix, const Vec (&x)[6][Batch][2],
                       Vec (&y)[Batch][2])
{
    const svbool_t pg = svptrue_b64();
    const Vec diagonal = svld1(pg, matrix + Row * lanes);
#pragma clang loop unroll(full)
    for (int b=0; b<Batch; ++b) {
        y[b][0] = svmul_x(pg, diagonal, x[Row][b][0]);
        y[b][1] = svmul_x(pg, diagonal, x[Row][b][1]);
    }
#pragma clang loop unroll(full)
    for (int col=0; col<6; ++col) {
        if (col == Row) continue;
        const bool upper = Row < col;
        const int idx = upper ? packed[Row][col] : packed[col][Row];
        const Vec ar = svld1(pg, matrix + idx * lanes);
        const Vec ai = svld1(pg, matrix + (idx+1) * lanes);
#pragma clang loop unroll(full)
        for (int b=0; b<Batch; ++b) {
            const Vec xr=x[col][b][0], xi=x[col][b][1];
            Vec yr=svmla_x(pg,y[b][0],ar,xr), yi=svmla_x(pg,y[b][1],ar,xi);
            y[b][0]=upper ? svmls_x(pg,yr,ai,xi) : svmla_x(pg,yr,ai,xi);
            y[b][1]=upper ? svmla_x(pg,yi,ai,xr) : svmls_x(pg,yi,ai,xr);
        }
    }
}

template<int Row, int Batch, bool Subtract>
inline void apply_row(const double* matrix, const Vec (&minus)[6][Batch][2],
                      const Vec (&plus)[6][Batch][2], double* const out[],
                      const double* const source[], int site)
{
    const svbool_t pg=svptrue_b64();
    Vec high[Batch][2], low[Batch][2];
    matvec_row<Row,Batch>(matrix,minus,high);
    matvec_row<Row,Batch>(matrix+36*lanes,plus,low);
    const int component=2*((Row/3)+4*(Row%3));
    const int offset=lanes*(spin_components*site+component);
#pragma clang loop unroll(full)
    for (int b=0; b<Batch; ++b) {
#pragma clang loop unroll(full)
        for (int ri=0; ri<2; ++ri) {
            Vec a=svsub_x(pg,low[b][ri],high[b][ri]);
            Vec c=svadd_x(pg,low[b][ri],high[b][ri]);
            if (Subtract) {
                a=svsub_x(pg,svld1(pg,source[b]+offset+lanes*ri),a);
                c=svsub_x(pg,svld1(pg,source[b]+offset+lanes*(ri+4)),c);
            }
            svst1(pg,out[b]+offset+lanes*ri,a);
            svst1(pg,out[b]+offset+lanes*(ri+4),c);
        }
    }
}

template<int Batch, bool Subtract>
inline void apply_site(int site,const double* clover,double* const out[],
                       const double* const in[],const double* const source[])
{
    const svbool_t pg=svptrue_b64();
    Vec minus[6][Batch][2], plus[6][Batch][2];
    const double* matrix=clover+site*clover_components*lanes;
    for (int row=0; row<6; ++row) {
        const int comp=2*(row/3+4*(row%3));
#pragma clang loop unroll(full)
        for (int b=0; b<Batch; ++b) {
            const double* x=in[b]+(site*spin_components+comp)*lanes;
            for (int ri=0; ri<2; ++ri) {
                const Vec a=svld1(pg,x+ri*lanes), c=svld1(pg,x+(ri+4)*lanes);
                minus[row][b][ri]=svsub_x(pg,c,a);
                plus[row][b][ri]=svadd_x(pg,c,a);
            }
        }
    }
    apply_row<0,Batch,Subtract>(matrix,minus,plus,out,source,site);
    apply_row<1,Batch,Subtract>(matrix,minus,plus,out,source,site);
    apply_row<2,Batch,Subtract>(matrix,minus,plus,out,source,site);
    apply_row<3,Batch,Subtract>(matrix,minus,plus,out,source,site);
    apply_row<4,Batch,Subtract>(matrix,minus,plus,out,source,site);
    apply_row<5,Batch,Subtract>(matrix,minus,plus,out,source,site);
}
template<int Batch, bool Subtract>
void apply_sites(int groups,const double* clover,double* const out[],
                 const double* const in[],const double* const source[])
{
#pragma omp for schedule(static)
    for (int site=0;site<groups;++site)
        apply_site<Batch,Subtract>(site,clover,out,in,source);
}
template<bool Subtract>
void tiled_single_sites(int groups,const double* clover,double* const out[],
                        const double* const in[],const double* const source[])
{
#pragma omp for schedule(static)
    for (int site=0;site<groups;++site) {
        for (int b=0;b<4;++b)
            apply_site<1,Subtract>(site,clover,out+b,in+b,Subtract ? source+b : nullptr);
    }
}

template<int Batch>
void apply(int groups,const double* clover,double* const out[],const double* const in[],
           const double* const source[])
{
    if (source) apply_sites<Batch,true>(groups,clover,out,in,source);
    else apply_sites<Batch,false>(groups,clover,out,in,nullptr);
}
}

namespace thesis_batch {
void tiled(int groups,const double* clover,double* const out[],
           const double* const in[],const double* const source[])
{
    if (source) tiled_single_sites<true>(groups,clover,out,in,source);
    else tiled_single_sites<false>(groups,clover,out,in,nullptr);
}
void inverse(int count,int groups,const double* clover,double* const out[],
             const double* const in[],const double* const source[])
{
    switch(count) {
        case 1: apply<1>(groups,clover,out,in,source); break;
        case 2: apply<2>(groups,clover,out,in,source); break;
        case 3: apply<3>(groups,clover,out,in,source); break;
        case 4: apply<4>(groups,clover,out,in,source); break;
        default: throw std::invalid_argument("batch count must be in [1,4]");
    }
}
}
#endif
