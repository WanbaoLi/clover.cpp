/*!
        @file    afield.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib_alt_QXS/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_QXS/inline/define_vlen.h"

#include "lib_alt_QXS/qxs_precision.h"


// template definition
#include "lib_alt_QXS/Field/afield-tmpl.h"


template<>
const std::string AField<real_t, QXS>::class_name = "AField<" QXS_PRECISION_NAME ", QXS>";


// explicit instanciation.
template class AField<real_t, QXS>;

//============================================================END=====
