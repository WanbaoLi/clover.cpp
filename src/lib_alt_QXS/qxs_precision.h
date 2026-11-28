#pragma once

#include "lib_alt_QXS/inline/define_vlen.h"

// Internal header for the sources compiled once per QXS precision.
#if QXS_SINGLE_PRECISION
#define VLEN VLENS
#define VLENX VLENXS
#define VLENY VLENYS
using real_t = float;
#define QXS_PRECISION_NAME "float"
#define QXS_WILSON_INTRINSICS "lib_alt_QXS/inline/vsimd_Wilson_SU3_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_float-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_float-inc.h"
#else
#define VLEN VLEND
#define VLENX VLENXD
#define VLENY VLENYD
using real_t = double;
#define QXS_PRECISION_NAME "double"
#define QXS_WILSON_INTRINSICS "lib_alt_QXS/inline/vsimd_Wilson_SU3_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_double-inc.h"
#include "lib_alt_QXS/inline/vsimd_common_double-inc.h"
#endif
