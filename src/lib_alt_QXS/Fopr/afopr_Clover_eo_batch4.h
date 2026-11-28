#pragma once
#include "configure.h"
#if defined(USE_ALT_QXS) && defined(USE_GROUP_SU3) && defined(__ARM_FEATURE_SVE_BITS) && __ARM_FEATURE_SVE_BITS == 512
#include "lib_alt_QXS/Fopr/afopr_Clover_eo.h"
#include "lib_alt_QXS/inline/define_vlen.h"
static_assert(VLEND == 8, "CloverEOBatch4 requires eight FP64 sites per pack");
#include <algorithm>
#include <stdexcept>

namespace thesis_batch {
void tiled(int groups,const double* clover,double* const out[],
           const double* const in[],const double* const source[]);
void inverse(int count,int groups,const double* clover,double* const out[],
             const double* const in[],const double* const source[]);

// Additive FP64 QXS operator. Construct and call in the same way as the existing
// AFopr_Clover_eo. The forward Schur methods are collective in the OpenMP team.
// Inputs, outputs and scratch are distinct checkerboard fields.
class CloverEOBatch4 : public AFopr_Clover_eo<AField<double,QXS>> {
public:
    using Field = AField<double,QXS>;
    using Base = AFopr_Clover_eo<Field>;
    explicit CloverEOBatch4(const Parameters& params) : Base(params) {}

    void local_inverse_batch(int count,Field* const out[],Field* const in[],
                             int parity,Field* const source[]=nullptr) {
        double* dst[4]; const double* src[4]; const double* sub[4];
        if (count<1 || count>4 || parity<0 || parity>1)
            throw std::invalid_argument("Invalid Clover batch size or parity");
        for (int b=0;b<count;++b) {
            dst[b]=out[b]->ptr(0); src[b]=in[b]->ptr(0);
            if (source) sub[b]=source[b]->ptr(0);
        }
        inverse(count,m_Nst2v,(parity ? m_To_inv : m_Te_inv).ptr(0),dst,src,
                source ? sub : nullptr);
    }

    void D_batch(int count,Field* const out[],Field* const in[],Field* const scratch[]) {
        if (count<1 || count>4) throw std::invalid_argument("Invalid Clover batch size");
        for (int b=0;b<count;++b) Meo(*scratch[b],*in[b],1);
        local_inverse_batch(count,out,scratch,1);
        for (int b=0;b<count;++b) Meo(*scratch[b],*out[b],0);
        local_inverse_batch(count,out,scratch,0,in);
    }

    void D_many(int nrhs,Field* const out[],Field* const in[],Field* const scratch[]) {
        if (nrhs<0) throw std::invalid_argument("Negative RHS count");
        for (int first=0; first<nrhs;) {
            const int count=std::min(nrhs-first,4);
            D_batch(count,out+first,in+first,scratch);
            first+=count;
        }
    }

    // Four independent width-one matvecs touch each packed site before moving on.
    // The last one to three columns use the checked width-specific inverse.
    void D_tiled(int nrhs,Field* const out[],Field* const in[],Field* const scratch[]) {
        if (nrhs<0) throw std::invalid_argument("Negative RHS count");
        int first=0;
        for (; nrhs-first>=4; first+=4)
            D_tiled_four(out+first,in+first,scratch);
        if (first<nrhs) D_batch(nrhs-first,out+first,in+first,scratch);
    }

private:
    void local_tiled(Field* const out[],Field* const in[],int parity,Field* const source[]=nullptr) {
        double* dst[4]; const double* src[4]; const double* sub[4];
        for (int b=0;b<4;++b) {
            dst[b]=out[b]->ptr(0); src[b]=in[b]->ptr(0);
            if (source) sub[b]=source[b]->ptr(0);
        }
        tiled(m_Nst2v,(parity ? m_To_inv : m_Te_inv).ptr(0),dst,src,source ? sub : nullptr);
    }
    void D_tiled_four(Field* const out[],Field* const in[],Field* const scratch[]) {
        for (int b=0;b<4;++b) Meo(*scratch[b],*in[b],1);
        local_tiled(out,scratch,1);
        for (int b=0;b<4;++b) Meo(*scratch[b],*out[b],0);
        local_tiled(out,scratch,0,in);
    }
};
}
#endif
