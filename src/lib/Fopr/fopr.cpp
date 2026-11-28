/*!
        @file    fopr.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2025-09-02 15:10:15 #$

        @version $LastChangedRevision: 2654 $
*/

#include "Fopr/fopr.h"
#include "Fopr/fopr_eo.h"
#include "Field/field.h"

#ifdef USE_FACTORY

#ifdef USE_FACTORY_AUTOREGISTER
#else
// setup factories for all subclasses

#include "fopr_Wilson.h"
#include "Org/fopr_Wilson_impl.h"
#include "Imp/fopr_Wilson_impl.h"
#include "fopr_Wilson_eo.h"
#include "Org/fopr_Wilson_eo_impl.h"
#include "Imp/fopr_Wilson_eo_impl.h"
#include "fopr_Clover.h"
#include "fopr_Clover_eo.h"
#include "fopr_Rational.h"
#include "fopr_Smeared.h"
#include "fopr_Smeared_eo.h"
#include "fopr_Chebyshev.h"

template<>
bool Fopr::init_factory()
{
  bool result = true;
  result &= Selector_Fopr_Wilson::register_factory();
  result &= Org::Fopr_Wilson::register_factory();
  result &= Imp::Fopr_Wilson::register_factory();
  result &= Selector_Fopr_Wilson_eo::register_factory();
  result &= Org::Fopr_Wilson_eo::register_factory();
  result &= Imp::Fopr_Wilson_eo::register_factory();
  result &= Fopr_Clover::register_factory();
  result &= Fopr_Clover_eo::register_factory();
  result &= Fopr_Rational::register_factory();
  result &= Fopr_Smeared::register_factory();
  result &= Fopr_Smeared_eo::register_factory();
  result &= Fopr_Chebyshev::register_factory();

  return result;
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
