/*!
      @file    bridgeQXS_Clover_coarse_float.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
      @date    $LastChangedDate: 2023-03-20 10:52:44 +0900 (2023年03月20日 (月)) $
      @version $LastChangedRevision: 2499 $
*/

#include "define_vlen.h"
#include "kernels/define_params.h"

#define  VLEN     VLENS
#define  VLENX    VLENXS
#define  VLENY    VLENYS

typedef float real_t;

//#include "define_index.h"
#include "vsimd_float-inc.h"
#include "vsimd_common_float-inc.h"
#include "vsimd_Wilson_SU3_float-inc.h"

#include "kernels/bridgeQXS_Clover_coarse.h"

#include "kernels/mult_Clover_coarse_parts_qxs-inc.h"

#include "kernels/mult_Clover_coarse_qxs-inc.h"
