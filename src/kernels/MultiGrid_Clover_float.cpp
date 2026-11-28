/*!

        @file    $Id: MultiGrid_Clover.cpp #$

        @brief   MultiGrid operation for Clover fermion (SIMD version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */

//====================================================================

#include "define_vlen.h"
#include "kernels/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#include "vsimd_float-inc.h"
#include "vsimd_common_float-inc.h"

#include "field/afield.h"
#include "field/afield-inc.h"
#include "field/afield_dd-inc.h"
#include "operator/afopr_Clover.h"
#include "operator/afopr_Clover_dd.h"
#include "field/aindex_coarse_lex.h"
#include "field/aindex_block_lex.h"

// template for MultiGrid_Clover
#include "solver/MultiGrid_Clover.h"
#include "solver/MultiGrid_Clover-tmpl.h"

typedef AField<float, QXS> AField_f;
//#define USE_IMPL_IN_TMPL

// specialization on single prec.


//====================================================================
template<>
const std::string MultiGrid_Clover<AField_f, AField_f>::class_name = "MultiGrid_Clover< AField<float,QXS>,  AField<float,QXS> >";
template class MultiGrid_Clover<AField_f, AField_f>;


//============================================================END=====
