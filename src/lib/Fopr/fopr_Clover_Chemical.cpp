/*!
        @file    fopr_Clover_Chemical.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$

        @version $LastChangedRevision: 2693 $
*/

#include "lib/Fopr/fopr_Clover_Chemical.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Clover_Chemical::register_factory();
}
#endif

#include "lib/Fopr/afopr_Clover_Chemical-tmpl.h"

template<>
const std::string Fopr_Clover_Chemical::class_name
                                         = "Fopr_Clover_Chemical";

template class AFopr_Clover_Chemical<Field>;

//============================================================END=====
