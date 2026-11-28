/*!
        @file    fopr.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-09-02 15:10:15 #$

        @version $LastChangedRevision: 2654 $
*/

#include "operator/fopr.h"
#include "operator/fopr_eo.h"
#include "field/field.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "operator/fopr_Clover.h"

template<>
bool Fopr::init_factory()
{
  return Fopr_Clover::register_factory();
}

#endif /* USE_FACTORY_AUTOREGISTER */

#endif /* USE_FACTORY */

// explicit instanciation.
template<>
const std::string AFopr<Field>::class_name = "Fopr";
template class AFopr<Field>;

template<>
const std::string AFopr_eo<Field>::class_name = "Fopr_eo";
template class AFopr_eo<Field>;
