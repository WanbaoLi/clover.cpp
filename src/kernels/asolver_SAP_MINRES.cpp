/*!
        @file    afopr_SAP_MINRES.cpp
        @brief   MINRES solver inside a SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2499 $
 */

//====================================================================
#include "solver/asolver_SAP_MINRES.h"
#include "runtime/threadManager.h"

#include "define_vlen.h"

#include "kernels/qxs_precision.h"

#include "field/aindex_block_lex.h"
#include "field/afield.h"
#include "field/afield-inc.h"
#include "field/afield_dd-inc.h"
#include "operator/afopr_Clover_dd.h"

//#define DEBUG_MINRES

template<typename AFIELD>
using AFopr_dd_t = AFopr_Clover_dd<AFIELD>;


#include "solver/asolver_SAP_MINRES-tmpl.h"

//====================================================================
template<>
const std::string ASolver_SAP_MINRES<AField<real_t, QXS> >::class_name
  = "ASolver_SAP_MINRES<AField<" QXS_PRECISION_NAME ",QXS> >";

template class ASolver_SAP_MINRES<AField<real_t, QXS> >;

//============================================================END=====
