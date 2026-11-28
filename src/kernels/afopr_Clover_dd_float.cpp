/*!
      @file    afopr_Clover_dd_float.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#include "operator/afopr_Clover_dd.h"

#include "runtime/threadManager.h"

#include "define_vlen.h"
#include "kernels/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

#define CHIRAL_ROTATION    // chiral rotation in clover term

#include "vsimd_float-inc.h"
#include "vsimd_common_float-inc.h"
#include "vsimd_Wilson_SU3_float-inc.h"

#include "field/aindex_lex.h"
#include "field/afield.h"
#include "field/afield-inc.h"
#include "field/afield_Gauge-inc.h"

#include "operator/mult_common_th-inc.h"
#include "operator/mult_Wilson_parts_qxs_org-inc.h"

#include "kernels/bridgeQXS_Wilson.h"
#include "kernels/bridgeQXS_Clover.h"


// template definition
#include "operator/afopr_Clover_dd-tmpl.h"

template<>
const std::string AFopr_Clover_dd<AField<float, QXS> >::class_name
  = "AFopr_Clover_dd<AField<float,QXS> >";

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init2 = AFopr<AField<float, QXS> >::Factory_params::Register(
    "Clover_dd", create_object_with_params);
}
#endif

// explicit instanciation
template class AFopr_Clover_dd<AField<float, QXS> >;

//============================================================END=====
