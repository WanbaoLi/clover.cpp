/*!
      @file    afopr_Clover.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Fopr/afopr_Clover.h"

#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt_QXS/inline/define_params.h"

#include "lib_alt_QXS/qxs_precision.h"

#define CHIRAL_ROTATION    // chiral rotation in clover term

#include QXS_WILSON_INTRINSICS

#include "lib_alt_QXS/Field/aindex_lex.h"
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Field/afield_Gauge-inc.h"

#include "lib_alt_QXS/Fopr/mult_common_th-inc.h"
#include "lib_alt_QXS/Fopr/mult_Wilson_parts_qxs_org-inc.h"

#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.h"
#include "lib_alt_QXS/BridgeQXS/bridgeQXS_Clover.h"


// template definition
#include "lib_alt_QXS/Fopr/afopr_Clover-tmpl.h"

template<>
const std::string AFopr_Clover<AField<real_t, QXS> >::class_name
  = "AFopr_Clover<AField<" QXS_PRECISION_NAME ",QXS> >";

// explicit instanciation
template class AFopr_Clover<AField<real_t, QXS> >;

//============================================================END=====
