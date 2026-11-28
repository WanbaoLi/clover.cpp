/*!

        @file    asolver_MG.cpp

        @brief   multigrid solver (QXS version)

        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate::  $

        @version $LastChangedRevision: 2499 $

 */
//====================================================================
#include "solver/asolver_MG.h"

#include "runtime/threadManager.h"
#include "runtime/randomNumberManager.h"
#include "field/afield.h"
#include "field/afield-inc.h"
#include "field/aindex_block_lex.h" // matching btw coasre/fine lattices

// clover specific
#include "operator/afopr_Clover.h"
#include "operator/afopr_Clover_dd.h"
#include "operator/afopr_Clover_coarse.h"
#include "solver/MultiGrid_Clover.h"


#include "solver/asolver_SAP_MINRES.h"
#include "solver/asolver_SAP.h"

typedef AField<float, QXS>    AField_f;
typedef AField<double, QXS>   AField_d;

// multigrid
using MultiGrid_t = MultiGrid_Clover<AField_f, AField_f>;

// operators
template<class Field> using FoprD_t = AFopr_Clover<Field>;
using FoprF_t      = AFopr_Clover_dd<AField_f>;
using FoprCoarse_t = AFopr_Clover_coarse<AField_f>;

// solver types
using CoarseSolver_t = ASolver_BiCGStab_Cmplx<AField_f>;
using Smoother_t = ASolver_SAP<AField_f>;

#include "solver/asolver_GCR.h"
#include "solver/asolver_MG-tmpl.h"

template<>
const std::string ASolver_MG<AField_d>::class_name = "ASolver_MG";

template class ASolver_MG<AField_d>;

template<>
const std::string ASolver_MG<AField_f>::class_name = "ASolver_MG";
template class ASolver_MG<AField_f>;
