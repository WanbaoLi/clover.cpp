/*!
        @file    afopr_SAP_MINRES.cpp
        @brief   MINRES solver inside a SAP solver (QXS version)
        @author  KANAMORI Issaku (kanamori)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate::  $
        @version $LastChangedRevision: 2499 $
 */

//====================================================================
#include "lib_alt/Solver/asolver_SAP_MINRES.h"
#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#include "lib_alt_QXS/qxs_precision.h"

#include "lib_alt_QXS/Field/aindex_block_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_dd-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"

//#define DEBUG_MINRES

template<typename AFIELD>
using AFopr_dd_t = AFopr_Clover_dd<AFIELD>;


#include "lib_alt/Solver/asolver_SAP_MINRES-tmpl.h"

//====================================================================
template<>
const std::string ASolver_SAP_MINRES<AField<real_t, QXS> >::class_name
  = "ASolver_SAP_MINRES<AField<" QXS_PRECISION_NAME ",QXS> >";

template class ASolver_SAP_MINRES<AField<real_t, QXS> >;

//============================================================END=====
