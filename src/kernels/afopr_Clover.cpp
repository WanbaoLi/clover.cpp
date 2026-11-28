/*!
      @file    afopr_Clover.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "operator/afopr_Clover.h"

#include "define_vlen.h"
#include "kernels/define_params.h"

#include "kernels/qxs_precision.h"

#define CHIRAL_ROTATION    // chiral rotation in clover term

#include QXS_WILSON_INTRINSICS

#include "field/aindex_lex.h"
#include "field/afield.h"
#include "field/afield-inc.h"
#include "field/afield_Gauge-inc.h"

#include "operator/mult_common_th-inc.h"
#include "operator/mult_Wilson_parts_qxs_org-inc.h"

#include "kernels/bridgeQXS_Wilson.h"
#include "kernels/bridgeQXS_Clover.h"


// template definition
#include "operator/afopr_Clover-tmpl.h"

template<>
const std::string AFopr_Clover<AField<real_t, QXS> >::class_name
  = "AFopr_Clover<AField<" QXS_PRECISION_NAME ",QXS> >";

// explicit instanciation
template class AFopr_Clover<AField<real_t, QXS> >;

//============================================================END=====
