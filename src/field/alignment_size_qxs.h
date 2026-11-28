#ifndef ALIGNMENT_SIZE_QXS_H
#define ALIGNMENT_SIZE_QXS_H

#include "field/alt_impl.h"
template<>
constexpr int alignment_size<QXS>() { return 256; }

#endif
